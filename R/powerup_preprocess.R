# R/powerup_preprocess.R

#' PowerUp preprocessing (platform contract)
#'
#' This function is intended to be called by the PowerUp Cloud Run worker.
#' It MUST read only local staged inputs and MUST write only under out_dir.
#'
#' @param job_id Character. Required job identifier for logging/manifest.
#' @param data_version Character. Dataset version tag (e.g. "v1").
#' @param response_set Character. "crispr" or "drug" (used for modelKey prefix).
#' @param seed Integer-like. Deterministic seed for hashing/splits.
#' @param local_input_matrix Path to local user matrix CSV (matrix.csv).
#' @param local_gene_expr Path to local gene_expression.csv.
#' @param local_response Path to local response CSV (CRISPR_response.csv or drug_response.csv staged as response.csv).
#' @param out_dir Path to write preprocess artifacts (e.g. /tmp/powerup/<job>/preprocess).
#' @param n_features Integer. Number of top-variance genes to select from user matrix (default 2000).
#' @param targets Optional character vector. Allowlisted perturbations to include (exact match).
#' @param sensitive_cutoff Numeric. Default 0.75.
#' @param resistant_cutoff Numeric. Default 0.25.
#'
#' @return Invisibly returns a list describing key outputs.
#' @export
powerup_preprocess <- function(
  job_id,
  data_version,
  response_set,
  seed,
  local_input_matrix,
  local_gene_expr,
  local_response,
  out_dir,
  n_features = 2000L,
  targets = NULL,
  sensitive_cutoff = 0.75,
  resistant_cutoff = 0.25
) {
  .pu_assert_scalar_chr(job_id, "job_id")
  .pu_assert_scalar_chr(data_version, "data_version")
  .pu_assert_scalar_chr(response_set, "response_set")
  .pu_assert_dir_writable(out_dir)

  seed_int <- .pu_as_int(seed, "seed")

  .pu_log(job_id, sprintf(
    "PREPROCESS start: data_version=%s response_set=%s seed=%d n_features=%d",
    data_version, response_set, seed_int, as.integer(n_features)
  ))

  # ---- Read inputs (strict) ----
  user_mat <- .pu_read_csv_strict(local_input_matrix, job_id, "local_input_matrix")
  gene_expr <- .pu_read_csv_strict(local_gene_expr, job_id, "local_gene_expr")
  resp <- .pu_read_csv_strict(local_response, job_id, "local_response")

  user_mat <- .pu_normalize_id_col(user_mat, "cell_line")
  gene_expr <- .pu_normalize_id_col(gene_expr, "cell_line")
  resp <- .pu_normalize_id_col(resp, "cell_line")

  # ---- Coerce numeric feature/response columns ----
  # Keep id col as character; coerce the rest to numeric (errors if impossible)
  user_mat <- .pu_coerce_numeric_cols(user_mat, id_col = "cell_line", job_id = job_id, label = "user_mat")
  gene_expr <- .pu_coerce_numeric_cols(gene_expr, id_col = "cell_line", job_id = job_id, label = "gene_expr")
  resp <- .pu_coerce_numeric_cols(resp, id_col = "cell_line", job_id = job_id, label = "response")

  # ---- Validate schema ----
  user_genes <- setdiff(colnames(user_mat), "cell_line")
  expr_genes <- setdiff(colnames(gene_expr), "cell_line")
  common_genes <- intersect(user_genes, expr_genes)

  if (length(common_genes) < 10) {
    .pu_fail(job_id, sprintf(
      "Too few overlapping genes between user matrix (%d) and gene_expression (%d): overlap=%d",
      length(user_genes), length(expr_genes), length(common_genes)
    ))
  }

  resp_cols <- setdiff(colnames(resp), "cell_line")
  if (length(resp_cols) < 1) {
    .pu_fail(job_id, "Response table has no perturbation columns (only cell_line).")
  }

  if (!is.null(targets)) {
    if (!is.character(targets) || length(targets) < 1) {
      .pu_fail(job_id, "targets must be a non-empty character vector when provided.")
    }
    unknown <- setdiff(targets, resp_cols)
    if (length(unknown) > 0) {
      .pu_fail(job_id, sprintf(
        "targets contains unknown perturbations: %s",
        paste(utils::head(unknown, 25), collapse = ", ")
      ))
    }
    resp_cols <- targets
  }

  # Train/test universe is intersection of gene_expr and resp rows
  all_ids <- intersect(gene_expr$cell_line, resp$cell_line)
  if (length(all_ids) < 20) {
    .pu_fail(job_id, sprintf("Too few overlapping cell_line IDs between gene_expression and response: n=%d", length(all_ids)))
  }
  all_ids <- sort(unique(all_ids))

  # Subset gene_expr/resp to the universe
  gene_expr_u <- gene_expr[match(all_ids, gene_expr$cell_line), , drop = FALSE]
  resp_u <- resp[match(all_ids, resp$cell_line), c("cell_line", resp_cols), drop = FALSE]

  # ---- Feature selection: top variable genes in user matrix (deterministic ties) ----
  n_features_int <- max(1L, as.integer(n_features))
  sel <- .pu_select_top_var_genes(user_mat, common_genes, n_features_int, job_id)
  selected_genes <- sel$genes
  actual_n_features <- length(selected_genes)

  # Build feature matrices
  features_all <- gene_expr_u[, c("cell_line", selected_genes), drop = FALSE]
  features_user <- user_mat[, c("cell_line", intersect(selected_genes, colnames(user_mat))), drop = FALSE]

  # Ensure user has all selected genes; if not, fail (keeps deterministic contract)
  missing_in_user <- setdiff(selected_genes, colnames(user_mat))
  if (length(missing_in_user) > 0) {
    .pu_fail(job_id, sprintf(
      "User matrix missing %d selected genes (unexpected after overlap filtering). Example: %s",
      length(missing_in_user), paste(utils::head(missing_in_user, 10), collapse = ", ")
    ))
  }

  # ---- Deterministic perturbation ordering + model keys ----
  perturbations <- sort(resp_cols)
  model_keys <- sprintf("%s_model_%05d", response_set, seq_along(perturbations))
  perturbations_df <- data.frame(
    modelKey = model_keys,
    perturbation = perturbations,
    stringsAsFactors = FALSE
  )

  # ---- Deterministic split selection (global test set) ----
  split <- .pu_build_split(
    job_id = job_id,
    seed_int = seed_int,
    ids = all_ids,
    response_df = resp_u,
    perturbations = perturbations,
    sensitive_cutoff = sensitive_cutoff,
    resistant_cutoff = resistant_cutoff
  )

  train_ids <- split$train_ids
  test_ids <- split$test_ids

  # ---- Materialize outcomes/features by split ----
  outcomes_train <- resp_u[match(train_ids, resp_u$cell_line), , drop = FALSE]
  outcomes_test  <- resp_u[match(test_ids,  resp_u$cell_line), , drop = FALSE]

  feat_train <- features_all[match(train_ids, features_all$cell_line), , drop = FALSE]
  feat_test  <- features_all[match(test_ids,  features_all$cell_line), , drop = FALSE]

  # Backward-compatible wide tables
  train_set <- cbind(outcomes_train, feat_train[, setdiff(colnames(feat_train), "cell_line"), drop = FALSE])
  test_set  <- cbind(outcomes_test,  feat_test[,  setdiff(colnames(feat_test),  "cell_line"), drop = FALSE])

  # ---- Write artifacts ----
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Required: perturbations.csv
  utils::write.csv(perturbations_df, file.path(out_dir, "perturbations.csv"), row.names = FALSE, quote = TRUE)

  # Better layout (recommended)
  utils::write.csv(feat_train, file.path(out_dir, "features_train.csv"), row.names = FALSE, quote = TRUE)
  utils::write.csv(feat_test,  file.path(out_dir, "features_test.csv"),  row.names = FALSE, quote = TRUE)
  utils::write.csv(features_user, file.path(out_dir, "features_user.csv"), row.names = FALSE, quote = TRUE)

  utils::write.csv(outcomes_train, file.path(out_dir, "outcomes_train.csv"), row.names = FALSE, quote = TRUE)
  utils::write.csv(outcomes_test,  file.path(out_dir, "outcomes_test.csv"),  row.names = FALSE, quote = TRUE)

  # Backward-compatible filenames (contract)
  utils::write.csv(train_set, file.path(out_dir, "train_set.csv"), row.names = FALSE, quote = TRUE)
  utils::write.csv(test_set,  file.path(out_dir, "test_set.csv"),  row.names = FALSE, quote = TRUE)
  utils::write.csv(features_user, file.path(out_dir, "user_samples.csv"), row.names = FALSE, quote = TRUE)

  # split.json
  .pu_write_json(split$split_json, file.path(out_dir, "split.json"))

  # manifest.json (minimum required fields + relative artifact paths)
  manifest <- list(
    schemaVersion = 1,
    jobId = job_id,
    createdAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
    seed = seed_int,
    dataVersion = data_version,
    responseSet = response_set,
    counts = list(
      totalModels = length(perturbations),
      nFeatures = actual_n_features,
      nTrain = length(train_ids),
      nTest = length(test_ids),
      nUser = nrow(features_user)
    ),
    artifacts = list(
      perturbationsCsv = "perturbations.csv",
      trainSetCsv = "train_set.csv",
      testSetCsv = "test_set.csv",
      userSamplesCsv = "user_samples.csv",
      featuresTrainCsv = "features_train.csv",
      featuresTestCsv = "features_test.csv",
      featuresUserCsv = "features_user.csv",
      outcomesTrainCsv = "outcomes_train.csv",
      outcomesTestCsv = "outcomes_test.csv",
      splitJson = "split.json"
    ),
    thresholds = list(
      sensitive = sensitive_cutoff,
      resistant = resistant_cutoff
    )
  )
  .pu_write_json(manifest, file.path(out_dir, "manifest.json"))

  .pu_log(job_id, sprintf(
    "PREPROCESS done: models=%d features=%d train=%d test=%d user=%d",
    length(perturbations), actual_n_features, length(train_ids), length(test_ids), nrow(features_user)
  ))

  invisible(list(
    out_dir = out_dir,
    total_models = length(perturbations),
    n_features = actual_n_features,
    n_train = length(train_ids),
    n_test = length(test_ids),
    n_user = nrow(features_user)
  ))
}

# ----------------------------
# Internal helpers (package-private)
# ----------------------------

.pu_log <- function(job_id, msg) {
  cat(sprintf("[powerup][jobId=%s] %s\n", job_id, msg))
}

.pu_fail <- function(job_id, msg) {
  stop(sprintf("[powerup][jobId=%s] %s", job_id, msg), call. = FALSE)
}

.pu_assert_scalar_chr <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || nchar(x) < 1) {
    stop(sprintf("Invalid %s: expected non-empty character scalar", name), call. = FALSE)
  }
}

.pu_as_int <- function(x, name) {
  if (is.integer(x) && length(x) == 1) return(x)
  if (is.numeric(x) && length(x) == 1) return(as.integer(x))
  if (is.character(x) && length(x) == 1 && grepl("^[0-9]+$", x)) return(as.integer(x))
  stop(sprintf("Invalid %s: expected integer-like scalar", name), call. = FALSE)
}

.pu_assert_dir_writable <- function(path) {
  if (!is.character(path) || length(path) != 1 || nchar(path) < 1) {
    stop("out_dir must be a non-empty string", call. = FALSE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  testfile <- file.path(path, ".write_test")
  ok <- tryCatch({
    cat("x", file = testfile)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) stop(sprintf("out_dir is not writable: %s", path), call. = FALSE)
  unlink(testfile, force = TRUE)
}

.pu_read_csv_strict <- function(path, job_id, label) {
  if (!file.exists(path)) .pu_fail(job_id, sprintf("Missing file for %s: %s", label, path))
  df <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) .pu_fail(job_id, sprintf("Failed reading %s CSV: %s", label, e$message))
  )
  if (!is.data.frame(df) || nrow(df) < 1 || ncol(df) < 2) {
    .pu_fail(job_id, sprintf("Invalid %s: expected >=1 row and >=2 columns: %s", label, path))
  }
  df
}

.pu_normalize_id_col <- function(df, expected = "cell_line") {
  # If expected exists, keep it; otherwise treat first column as id and rename
  if (!(expected %in% colnames(df))) {
    colnames(df)[1] <- expected
  }
  df[[expected]] <- as.character(df[[expected]])
  df
}

.pu_coerce_numeric_cols <- function(df, id_col, job_id, label) {
  cols <- setdiff(colnames(df), id_col)
  if (length(cols) < 1) .pu_fail(job_id, sprintf("%s has no numeric columns besides %s", label, id_col))
  for (cn in cols) {
    v <- df[[cn]]
    if (is.character(v)) {
      suppressWarnings(vn <- as.numeric(v))
      if (all(is.na(vn)) && any(nzchar(v))) {
        .pu_fail(job_id, sprintf("%s column '%s' is non-numeric and cannot be coerced", label, cn))
      }
      df[[cn]] <- vn
    } else if (!is.numeric(v)) {
      suppressWarnings(vn <- as.numeric(v))
      df[[cn]] <- vn
    }
  }
  df
}

.pu_select_top_var_genes <- function(user_mat, gene_candidates, n_features, job_id) {
  # user_mat: data.frame with cell_line + genes
  genes <- intersect(gene_candidates, colnames(user_mat))
  if (length(genes) < 10) .pu_fail(job_id, "Too few gene candidates for variance selection.")
  m <- as.matrix(user_mat[, genes, drop = FALSE])
  # variance per column (fast + deterministic)
  vars <- apply(m, 2, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- -Inf

  ord <- order(-vars, names(vars))  # tie-break by gene name asc
  genes_sorted <- names(vars)[ord]
  keep <- utils::head(genes_sorted, n_features)
  keep <- keep[is.finite(vars[keep])]

  if (length(keep) < 10) .pu_fail(job_id, "Variance selection yielded too few usable genes after filtering.")
  list(genes = keep, variances = vars[keep])
}

.pu_hash_to_index <- function(seed_int, perturbation, klass, n) {
  # Deterministic, platform-independent hash -> 1..n
  # Uses digest::digest; add digest to Imports.
  key <- paste(seed_int, perturbation, klass, sep = "|")
  h <- digest::digest(key, algo = "xxhash32", serialize = FALSE)
  # take first 8 hex chars -> integer
  x <- strtoi(substr(h, 1, 8), base = 16)
  (x %% n) + 1L
}

.pu_build_split <- function(job_id, seed_int, ids, response_df, perturbations,
                            sensitive_cutoff, resistant_cutoff) {
  # response_df includes cell_line + perturbations
  # ids: sorted vector of cell_line
  # Build global test set by selecting one sensitive + one resistant per perturbation
  test_set <- character(0)
  decisions <- vector("list", length(perturbations))
  names(decisions) <- perturbations

  # Pre-map row index by id for quick subsetting
  idx_map <- setNames(seq_len(nrow(response_df)), response_df$cell_line)

  for (p in perturbations) {
    y <- response_df[[p]]
    if (is.null(y)) next

    # Candidates among ids only (response_df is already aligned, but keep explicit)
    # Sensitive: >= cutoff; Resistant: <= cutoff
    sens_ids <- response_df$cell_line[which(!is.na(y) & y >= sensitive_cutoff)]
    res_ids  <- response_df$cell_line[which(!is.na(y) & y <= resistant_cutoff)]

    sens_ids <- intersect(sort(unique(sens_ids)), ids)
    res_ids  <- intersect(sort(unique(res_ids)), ids)

    pick_one <- function(cands, klass) {
      if (length(cands) == 0) return(NA_character_)
      start <- .pu_hash_to_index(seed_int, p, klass, length(cands))
      # collision-resistant deterministic stepping
      for (k in seq_len(length(cands))) {
        j <- ((start - 1L + (k - 1L)) %% length(cands)) + 1L
        choice <- cands[j]
        if (!(choice %in% test_set)) return(choice)
      }
      # If all collide, still return deterministic first
      cands[start]
    }

    sens_pick <- pick_one(sens_ids, "sensitive")
    if (!is.na(sens_pick)) test_set <- c(test_set, sens_pick)

    res_pick <- pick_one(res_ids, "resistant")
    if (!is.na(res_pick)) test_set <- c(test_set, res_pick)

    decisions[[p]] <- list(
      nSensitive = length(sens_ids),
      nResistant = length(res_ids),
      pickedSensitive = sens_pick,
      pickedResistant = res_pick
    )
  }

  test_set <- sort(unique(test_set))
  train_set <- setdiff(ids, test_set)

  if (length(test_set) < 2) {
    .pu_fail(job_id, "Split produced too small of a test set; check response thresholds / data.")
  }
  if (length(train_set) < 10) {
    .pu_fail(job_id, "Split produced too small of a train set.")
  }

  split_json <- list(
    schemaVersion = 1,
    seed = seed_int,
    sensitiveCutoff = sensitive_cutoff,
    resistantCutoff = resistant_cutoff,
    trainIds = train_set,
    testIds = test_set,
    perPerturbation = decisions
  )

  list(
    train_ids = train_set,
    test_ids = test_set,
    split_json = split_json
  )
}

.pu_write_json <- function(obj, path) {
  jsonlite::write_json(obj, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}