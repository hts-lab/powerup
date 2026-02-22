# R/powerup_platform_api.R

#' POWERUP platform API (Cloud Run + standalone)
#'
#' These functions are called by the POWERUP worker container entrypoint.
#' They MUST remain deterministic given the same inputs + seed.
#'
#' @name powerup_platform_api
#' @keywords internal
NULL

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tibble)
  library(glue)
  library(digest)
  library(data.table)
})

# ---- small utilities ----

powerup_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

powerup_write_json <- function(path, obj) {
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE), con = path)
}

# ---- internal helpers (package-private) ----

.pu_assert_file_exists <- function(path, label) {
  if (!is.character(path) || length(path) != 1 || nchar(path) < 1) {
    stop(glue("{label} must be a non-empty string"))
  }
  if (!file.exists(path)) {
    stop(glue("{label} does not exist: {path}"))
  }
}

.pu_normalize_id_col <- function(df, expected = "cell_line") {
  if (!(expected %in% colnames(df))) {
    # If the file doesn't have cell_line explicitly, treat first column as id
    colnames(df)[1] <- expected
  }
  df[[expected]] <- as.character(df[[expected]])
  df
}

.pu_coerce_numeric_cols <- function(df, id_col, label) {
  cols <- setdiff(colnames(df), id_col)
  if (length(cols) < 1) stop(glue("{label} has no columns besides {id_col}"))

  for (cn in cols) {
    v <- df[[cn]]
    if (is.character(v)) {
      suppressWarnings(vn <- as.numeric(v))
      # If there were non-empty strings but coercion yields all NA -> hard fail
      if (all(is.na(vn)) && any(nzchar(v))) {
        stop(glue("{label} column '{cn}' is non-numeric and cannot be coerced"))
      }
      df[[cn]] <- vn
    } else if (!is.numeric(v)) {
      suppressWarnings(df[[cn]] <- as.numeric(v))
    }
  }
  df
}

.pu_select_top_var_genes <- function(user_matrix, genes, n_features) {
  genes <- intersect(genes, colnames(user_matrix))
  m <- as.matrix(user_matrix[, genes, drop = FALSE])

  vars <- apply(m, 2, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- -Inf

  ord <- order(-vars, names(vars))   # deterministic tiebreak by gene name
  keep <- names(vars)[ord]
  keep <- keep[seq_len(min(length(keep), n_features))]
  keep <- keep[is.finite(vars[keep])]
  keep
}

.pu_hash_to_index <- function(seed_int, perturbation, klass, n) {
  key <- paste(seed_int, perturbation, klass, sep = "|")
  h <- digest::digest(key, algo = "xxhash32", serialize = FALSE)

  # Use 7 hex chars to avoid overflow/NA from strtoi()
  x <- strtoi(substr(h, 1, 7), base = 16)
  if (is.na(x) || n <= 0) stop("hash_to_index produced NA or invalid n")

  (x %% n) + 1L
}

.pu_build_split <- function(
  seed_int,
  response_u,
  perturbations,
  sensitive_cutoff,
  resistant_cutoff,
  n_per_class = 1L,      # kept for signature compatibility (we target >=1 each class under budget)
  max_test_frac = 0.25   # hard cap for test set fraction as % of total cell lines
) {
  # -----------------------------
  # Setup + sanity
  # -----------------------------
  perturbations <- as.character(perturbations)
  perturbations <- perturbations[perturbations %in% colnames(response_u)]
  perturbations <- sort(unique(perturbations))

  all_ids <- sort(unique(as.character(response_u$cell_line)))
  n_all <- length(all_ids)
  if (n_all < 20) stop("Too few cell lines to split")

  budget <- max(2L, min(n_all - 10L, as.integer(ceiling(max_test_frac * n_all))))
  if (budget <= 1L) stop("Invalid test budget computed")

  # Needs: for each perturbation, we want >=1 sensitive and >=1 resistant in test
  need_sens <- rep(TRUE, length(perturbations))
  need_res  <- rep(TRUE, length(perturbations))
  names(need_sens) <- perturbations
  names(need_res)  <- perturbations

  # Deterministic tiebreak hash for a cell line
  hash_rank <- function(cell_line) {
    h <- digest::digest(paste(seed_int, cell_line, sep = "|"), algo = "xxhash32", serialize = FALSE)
    x <- strtoi(substr(h, 1, 7), base = 16)
    if (is.na(x)) 0L else x
  }

  # -----------------------------
  # Precompute coverage sets per cell line
  #   sens_cov[[i]] = indices of perturbations where cell_line i is sensitive
  #   res_cov[[i]]  = indices of perturbations where cell_line i is resistant
  # -----------------------------
  # Map ids to row index (response_u already filtered to overlap universe upstream)
  id_vec <- as.character(response_u$cell_line)

  # We'll build coverage via one pass over perturbations (column-wise)
  n_p <- length(perturbations)
  sens_cov <- vector("list", n_all)
  res_cov  <- vector("list", n_all)
  names(sens_cov) <- all_ids
  names(res_cov)  <- all_ids
  for (i in seq_len(n_all)) {
    sens_cov[[i]] <- integer(0)
    res_cov[[i]]  <- integer(0)
  }

  # To avoid repeated matching, precompute row indices for each all_id
  # (response_u may already be in same order as all_ids; but we do a safe match)
  row_of_id <- match(all_ids, id_vec)
  if (any(is.na(row_of_id))) {
    # Shouldn't happen since all_ids derived from response_u, but keep robust
    row_of_id <- which(id_vec %in% all_ids)
  }

  # Column-wise scan
  for (j in seq_len(n_p)) {
    p <- perturbations[[j]]
    y <- response_u[[p]]

    # Pull values for all_ids in stable order
    v <- y[row_of_id]

    sens_mask <- !is.na(v) & (v >= sensitive_cutoff)
    res_mask  <- !is.na(v) & (v <= resistant_cutoff)

    if (any(sens_mask)) {
      ids <- which(sens_mask)
      for (ii in ids) sens_cov[[ii]] <- c(sens_cov[[ii]], j)
    }
    if (any(res_mask)) {
      ids <- which(res_mask)
      for (ii in ids) res_cov[[ii]] <- c(res_cov[[ii]], j)
    }
  }

  # -----------------------------
  # Greedy selection under budget
  # -----------------------------
  selected <- rep(FALSE, n_all)
  test_set <- character(0)

  # Track unmet needs as logical vectors for fast sum()
  need_sens_vec <- rep(TRUE, n_p)
  need_res_vec  <- rep(TRUE, n_p)

  # Helper: compute marginal gain of candidate i
  marginal_gain <- function(i) {
    if (selected[[i]]) return(-1L)

    s_idx <- sens_cov[[i]]
    r_idx <- res_cov[[i]]

    gain <- 0L
    if (length(s_idx) > 0) gain <- gain + sum(need_sens_vec[s_idx])
    if (length(r_idx) > 0) gain <- gain + sum(need_res_vec[r_idx])

    gain
  }

  # For determinism and speed, we’ll iterate and recompute gains each round.
  # n_all ~1121 and budget ~280 => ~300k gain calcs; fine.
  for (step in seq_len(budget)) {
    # Stop early if all needs met
    if (!any(need_sens_vec) && !any(need_res_vec)) break

    best_i <- NA_integer_
    best_gain <- -1L
    best_hash <- Inf

    for (i in seq_len(n_all)) {
      if (selected[[i]]) next
      g <- marginal_gain(i)
      if (g < best_gain) next

      # tie-break deterministically
      h <- hash_rank(all_ids[[i]])
      if (g > best_gain || (g == best_gain && h < best_hash)) {
        best_gain <- g
        best_i <- i
        best_hash <- h
      }
    }

    # If no candidate adds anything, stop
    if (is.na(best_i) || best_gain <= 0L) break

    # Select it
    selected[[best_i]] <- TRUE
    test_set <- c(test_set, all_ids[[best_i]])

    # Update unmet needs
    s_idx <- sens_cov[[best_i]]
    r_idx <- res_cov[[best_i]]
    if (length(s_idx) > 0) need_sens_vec[s_idx] <- FALSE
    if (length(r_idx) > 0) need_res_vec[r_idx]  <- FALSE
  }

  test_set <- sort(unique(test_set))
  train_set <- setdiff(all_ids, test_set)

  n_test  <- length(test_set)
  n_train <- length(train_set)

  # -----------------------------
  # Build decisions (compact!)
  # For each perturbation, store ONE representative pickedSensitive/pickedResistant
  # from within the final test_set (deterministic: lexicographically smallest id).
  # -----------------------------
  test_idx <- match(test_set, all_ids)
  decisions <- vector("list", n_p)
  names(decisions) <- perturbations

  # Precompute membership for fast lookups
  test_member <- rep(FALSE, n_all)
  test_member[test_idx] <- TRUE

  for (j in seq_len(n_p)) {
    p <- perturbations[[j]]
    y <- response_u[[p]][row_of_id]

    sens_mask <- test_member & !is.na(y) & (y >= sensitive_cutoff)
    res_mask  <- test_member & !is.na(y) & (y <= resistant_cutoff)

    sens_ids <- all_ids[which(sens_mask)]
    res_ids  <- all_ids[which(res_mask)]

    picked_s <- if (length(sens_ids) > 0) sort(sens_ids)[[1]] else NA_character_
    picked_r <- if (length(res_ids) > 0) sort(res_ids)[[1]] else NA_character_

    decisions[[j]] <- list(
      nSensitive = sum(!is.na(y) & (y >= sensitive_cutoff)),
      nResistant = sum(!is.na(y) & (y <= resistant_cutoff)),
      pickedSensitive = picked_s,
      pickedResistant = picked_r
    )
  }

  # Coverage stats
  has_s <- vapply(decisions, function(d) !is.na(d$pickedSensitive), logical(1))
  has_r <- vapply(decisions, function(d) !is.na(d$pickedResistant), logical(1))
  has_both <- has_s & has_r

  cat(sprintf(
    paste0(
      "[powerup] split summary (budgeted greedy): ",
      "nPerturbations=%d nAllIds=%d budget=%d nTest=%d (%.1f%%) nTrain=%d | ",
      "coveredSensitive=%d (%.1f%%) coveredResistant=%d (%.1f%%) coveredBoth=%d (%.1f%%)\n"
    ),
    n_p, n_all, budget, n_test, 100.0 * n_test / n_all, n_train,
    sum(has_s), 100.0 * mean(has_s),
    sum(has_r), 100.0 * mean(has_r),
    sum(has_both), 100.0 * mean(has_both)
  ))

  # Hard constraints
  if (n_test < 2) {
    stop(sprintf("Split produced too small test set | nAllIds=%d nTest=%d nTrain=%d nPerturbations=%d", n_all, n_test, n_train, n_p))
  }
  if (n_train < 10) {
    stop(sprintf(
      "Split produced too small train set | nAllIds=%d nTest=%d nTrain=%d nPerturbations=%d | firstTestIds=%s",
      n_all, n_test, n_train, n_p, paste(head(test_set, 10), collapse = ",")
    ))
  }
  if (n_test > budget) {
    stop(sprintf(
      "Split violated test budget | nAllIds=%d budget=%d nTest=%d (%.1f%%)",
      n_all, budget, n_test, 100.0 * n_test / n_all
    ))
  }

  split_json <- list(
    schemaVersion = 2,
    seed = seed_int,
    sensitiveCutoff = sensitive_cutoff,
    resistantCutoff = resistant_cutoff,
    maxTestFrac = max_test_frac,
    budget = budget,
    trainIds = train_set,
    testIds = test_set,
    coverage = list(
      coveredSensitive = sum(has_s),
      coveredResistant = sum(has_r),
      coveredBoth = sum(has_both),
      totalPerturbations = n_p
    ),
    perPerturbation = decisions
  )

  list(
    train_ids = train_set,
    test_ids = test_set,
    split_json = split_json
  )
}

# ---- CONTRACT 1: preprocess ----
# Deterministic preprocessing owned by the R package.
# Reads local staged inputs, writes only under out_preprocess_dir.
#' @export
powerup_preprocess <- function(
  gene_expression_path,
  response_path,
  matrix_path,
  out_preprocess_dir,
  data_version,
  response_set,
  seed,
  job_id,
  n_features = 2000L,
  targets = NULL,
  sensitive_cutoff = 0.75,
  resistant_cutoff = 0.25
) {
  powerup_dir_create(out_preprocess_dir)

  seed_int <- as.integer(seed)
  if (is.na(seed_int)) stop("seed must be integer-like")

  message(glue("[powerup][jobId={job_id}] PREPROCESS start data_version={data_version} response_set={response_set} seed={seed_int}"))

  # ---- Read inputs (strict) ----
  .pu_assert_file_exists(gene_expression_path, "gene_expression_path")
  .pu_assert_file_exists(response_path, "response_path")
  .pu_assert_file_exists(matrix_path, "matrix_path")

  gene_expression <- readr::read_csv(gene_expression_path, show_col_types = FALSE, progress = FALSE)
  response_df     <- readr::read_csv(response_path, show_col_types = FALSE, progress = FALSE)
  user_matrix     <- readr::read_csv(matrix_path, show_col_types = FALSE, progress = FALSE)

  # Normalize ID column name to cell_line (if first col not named)
  gene_expression <- .pu_normalize_id_col(gene_expression, "cell_line")
  response_df     <- .pu_normalize_id_col(response_df, "cell_line")
  user_matrix     <- .pu_normalize_id_col(user_matrix, "cell_line")

  # Coerce numerics (fail fast if non-numeric)
  gene_expression <- .pu_coerce_numeric_cols(gene_expression, id_col = "cell_line", label = "gene_expression")
  response_df     <- .pu_coerce_numeric_cols(response_df,     id_col = "cell_line", label = "response_df")
  user_matrix     <- .pu_coerce_numeric_cols(user_matrix,     id_col = "cell_line", label = "user_matrix")

  # ---- Validate overlaps ----
  user_genes <- setdiff(colnames(user_matrix), "cell_line")
  expr_genes <- setdiff(colnames(gene_expression), "cell_line")
  common_genes <- intersect(user_genes, expr_genes)
  if (length(common_genes) < 10) {
    stop(glue("[powerup][jobId={job_id}] Too few overlapping genes between user matrix and gene_expression: overlap={length(common_genes)}"))
  }

  perturbations <- setdiff(colnames(response_df), "cell_line")
  if (length(perturbations) < 1) stop(glue("[powerup][jobId={job_id}] response_df has no perturbation columns"))

  if (!is.null(targets)) {
    if (!is.character(targets) || length(targets) < 1) stop(glue("[powerup][jobId={job_id}] targets must be non-empty character vector"))
    unknown <- setdiff(targets, perturbations)
    if (length(unknown) > 0) stop(glue("[powerup][jobId={job_id}] targets contains unknown perturbations: {paste(head(unknown,25), collapse=', ')}"))
    perturbations <- targets
  }

  # Universe = cell_line overlap between gene_expression and response
  all_ids <- intersect(gene_expression$cell_line, response_df$cell_line)
  all_ids <- sort(unique(all_ids))
  if (length(all_ids) < 20) stop(glue("[powerup][jobId={job_id}] Too few overlapping cell_line IDs between gene_expression and response: n={length(all_ids)}"))

  gene_expression_u <- gene_expression %>% dplyr::filter(.data$cell_line %in% all_ids)
  response_u        <- response_df %>%
    dplyr::filter(.data$cell_line %in% all_ids) %>%
    dplyr::select(.data$cell_line, dplyr::all_of(perturbations))

  # ---- Feature selection: top variable genes on user matrix (deterministic ties) ----
  n_features_int <- max(1L, as.integer(n_features))
  selected_genes <- .pu_select_top_var_genes(user_matrix, common_genes, n_features_int)
  if (length(selected_genes) < 10) stop(glue("[powerup][jobId={job_id}] Feature selection produced too few genes: n={length(selected_genes)}"))

  # Build feature tables
  features_all <- gene_expression_u %>% dplyr::select(.data$cell_line, dplyr::all_of(selected_genes))
  features_user <- user_matrix %>% dplyr::select(.data$cell_line, dplyr::all_of(selected_genes))

  # ---- Deterministic perturbations table + model keys ----
  perturbations <- sort(perturbations)
  model_keys <- sprintf("%s_model_%05d", response_set, seq_along(perturbations))
  perturbations_tbl <- tibble::tibble(
    modelKey = model_keys,
    perturbation = perturbations
  )

  # ---- Deterministic split (global test set) ----
  split <- .pu_build_split(
    seed_int = seed_int,
    response_u = response_u,
    perturbations = perturbations,
    sensitive_cutoff = sensitive_cutoff,
    resistant_cutoff = resistant_cutoff,
    n_per_class = 1L
  )

  train_ids <- split$train_ids
  test_ids  <- split$test_ids

  # Materialize outcomes/features by split
  outcomes_train <- response_u %>% dplyr::filter(.data$cell_line %in% train_ids)
  outcomes_test  <- response_u %>% dplyr::filter(.data$cell_line %in% test_ids)

  features_train <- features_all %>% dplyr::filter(.data$cell_line %in% train_ids)
  features_test  <- features_all %>% dplyr::filter(.data$cell_line %in% test_ids)

  # ---- Write artifacts (stable filenames) ----
  readr::write_csv(perturbations_tbl, file.path(out_preprocess_dir, "perturbations.csv"))

  # Better split layout ONLY (no legacy wide tables)
  readr::write_csv(features_train, file.path(out_preprocess_dir, "features_train.csv"))
  readr::write_csv(features_test,  file.path(out_preprocess_dir, "features_test.csv"))
  readr::write_csv(features_user,  file.path(out_preprocess_dir, "features_user.csv"))
  readr::write_csv(outcomes_train, file.path(out_preprocess_dir, "outcomes_train.csv"))
  readr::write_csv(outcomes_test,  file.path(out_preprocess_dir, "outcomes_test.csv"))

  # split.json
  powerup_write_json(file.path(out_preprocess_dir, "split.json"), split$split_json)

  # manifest.json (minimum required)
  manifest <- list(
    schemaVersion = 1,
    jobId = job_id,
    createdAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
    seed = seed_int,
    dataVersion = data_version,
    responseSet = response_set,
    counts = list(
      totalModels = length(perturbations),
      nFeatures = length(selected_genes),
      nTrain = length(train_ids),
      nTest = length(test_ids),
      nUser = nrow(features_user)
    ),
    artifacts = list(
      perturbationsCsv = "perturbations.csv",
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
  powerup_write_json(file.path(out_preprocess_dir, "manifest.json"), manifest)

  message(glue("[powerup][jobId={job_id}] PREPROCESS done models={length(perturbations)} features={length(selected_genes)} train={length(train_ids)} test={length(test_ids)} user={nrow(features_user)}"))
  invisible(TRUE)
}

# ---- CONTRACT 2: train models ----
# This function receives a list of model keys and trains them.
#' @export
powerup_train_models <- function(
  train_set_path,
  test_set_path,
  user_samples_path,
  perturbations_path,
  model_keys,
  out_models_dir,
  seed,
  response_set,
  job_id,

  # Defaults for make_xgb_model()
  response_cutoff = 0.75,
  decreasing = FALSE,
  weight_cap = 0,
  nfolds = 3,
  nrepeats = 3,
  nrounds = 100,
  max_depth = 3,
  f_subsample = 1,
  min_score = 0.05,
  skip_eval = FALSE,
  shuffle = FALSE,
  n_threads = 4,
  xgb_params = NULL,
  cor_data = NULL,
  cor_n_features = 1000,
  use_gpu = FALSE,
  gpu_id = 0
) {
  powerup_dir_create(out_models_dir)
  set.seed(as.integer(seed))

  # ---- Validate model_keys ----
  if (is.null(model_keys) || length(model_keys) == 0) {
    stop(glue("[powerup][jobId={job_id}] model_keys is empty"))
  }
  if (!is.character(model_keys)) {
    stop(glue("[powerup][jobId={job_id}] model_keys must be a character vector"))
  }
  if (anyNA(model_keys) || any(trimws(model_keys) == "")) {
    bad <- model_keys[is.na(model_keys) | trimws(model_keys) == ""]
    stop(glue("[powerup][jobId={job_id}] model_keys contains NA/empty values: {paste(bad, collapse=', ')}"))
  }
  if (any(duplicated(model_keys))) {
    dups <- unique(model_keys[duplicated(model_keys)])
    stop(glue("[powerup][jobId={job_id}] model_keys contains duplicates: {paste(dups, collapse=', ')}"))
  }

  message(glue("[powerup][jobId={job_id}] Received {length(model_keys)} model_keys (sample: {paste(head(model_keys, 10), collapse=', ')})"))


  # ---- perturbations table is the mapping source of truth ----
  pert_tbl <- readr::read_csv(perturbations_path, show_col_types = FALSE, progress = FALSE)
  message(glue("[powerup][jobId={job_id}] perturbations_path={perturbations_path}"))
  message(glue("[powerup][jobId={job_id}] perturbations.csv rows={nrow(pert_tbl)} cols={ncol(pert_tbl)}"))
  message(glue("[powerup][jobId={job_id}] perturbations.csv colnames={paste(names(pert_tbl), collapse=', ')}"))
  if (!("modelKey" %in% names(pert_tbl))) stop("perturbations.csv must contain column: modelKey")
  if (!("perturbation" %in% names(pert_tbl))) {
    stop("perturbations.csv must contain column: perturbation (maps modelKey -> dataset outcome column)")
  }


  key_map <- pert_tbl %>%
    filter(.data$modelKey %in% model_keys) %>%
    select(.data$modelKey, .data$perturbation)

  if (nrow(key_map) != length(model_keys)) {
    missing <- setdiff(model_keys, key_map$modelKey)
    stop(glue("Some modelKeys missing in perturbations.csv: {paste(missing, collapse=', ')}"))
  }

  # Check perturbation values look sane
  if (any(is.na(key_map$perturbation) | trimws(key_map$perturbation) == "")) {
    bad <- key_map$modelKey[is.na(key_map$perturbation) | trimws(key_map$perturbation) == ""]
    stop(glue("[powerup][jobId={job_id}] Some mapped perturbations are NA/empty for modelKeys: {paste(bad, collapse=', ')}"))
  }

  # Write a shard-level debug artifact so you can inspect exactly what was trained
  debug_map_path <- file.path(out_models_dir, "_debug_key_map.csv")
  readr::write_csv(key_map, debug_map_path)
  message(glue("[powerup][jobId={job_id}] Wrote debug mapping: {debug_map_path}"))
  message(glue("[powerup][jobId={job_id}] key_map sample:\n{paste(capture.output(print(utils::head(key_map, 10))), collapse='\n')}"))



  # ---- Split layout ONLY ----
  # NOTE: train_set_path/test_set_path/user_samples_path are retained for signature compatibility,
  # but we interpret them as anchors to find split-layout artifacts in the same directory.
  features_train_path <- file.path(dirname(train_set_path), "features_train.csv")
  features_test_path  <- file.path(dirname(test_set_path),  "features_test.csv")
  features_user_path  <- file.path(dirname(user_samples_path), "features_user.csv")
  outcomes_train_path <- file.path(dirname(train_set_path), "outcomes_train.csv")
  outcomes_test_path  <- file.path(dirname(test_set_path),  "outcomes_test.csv")

  .pu_assert_file_exists(features_train_path, "features_train.csv")
  .pu_assert_file_exists(features_test_path, "features_test.csv")
  .pu_assert_file_exists(features_user_path, "features_user.csv")
  .pu_assert_file_exists(outcomes_train_path, "outcomes_train.csv")
  .pu_assert_file_exists(outcomes_test_path, "outcomes_test.csv")

  message(glue("[powerup][jobId={job_id}] Using split layout ONLY (features_* + outcomes_* with column-pruning)"))

  # Features (shared across all models)
  feat_train <- readr::read_csv(features_train_path, show_col_types = FALSE, progress = FALSE) %>%
    tibble::column_to_rownames("cell_line")
  feat_test <- readr::read_csv(features_test_path, show_col_types = FALSE, progress = FALSE) %>%
    tibble::column_to_rownames("cell_line")
  user_df <- readr::read_csv(features_user_path, show_col_types = FALSE, progress = FALSE) %>%
    tibble::column_to_rownames("cell_line")

  # Outcomes: read only the needed columns for this shard (cell_line + perturbations in shard)
  needed_outcomes <- unique(key_map$perturbation)

  out_train_tbl <- data.table::fread(
    outcomes_train_path,
    select = c("cell_line", needed_outcomes),
    data.table = FALSE
  )
  out_test_tbl <- data.table::fread(
    outcomes_test_path,
    select = c("cell_line", needed_outcomes),
    data.table = FALSE
  )

  out_train <- out_train_tbl %>% tibble::column_to_rownames("cell_line")
  out_test  <- out_test_tbl  %>% tibble::column_to_rownames("cell_line")

  # We train each model and write artifacts per modelKey
  total <- length(model_keys)

  for (i in seq_along(model_keys)) {
    mk <- model_keys[[i]]
    perturbation <- key_map$perturbation[key_map$modelKey == mk][[1]]

    model_out_dir <- file.path(out_models_dir, mk)
    powerup_dir_create(model_out_dir)

    # Per-model guardrails: never let a single failure kill the shard.
    # We still write model artifacts with error details so the platform can proceed.
    tryCatch({

      message(glue("[powerup][jobId={job_id}] TRAIN modelKey={mk} perturbation={perturbation} ({i} of {total}) start"))

      if (!(perturbation %in% colnames(out_train))) {
        stop(glue("[powerup][jobId={job_id}] Missing outcome column '{perturbation}' in outcomes_train.csv"))
      }
      if (!(perturbation %in% colnames(out_test))) {
        stop(glue("[powerup][jobId={job_id}] Missing outcome column '{perturbation}' in outcomes_test.csv"))
      }

      # Build per-model dataset: outcome column + all features
      train_df <- cbind(
        setNames(as.data.frame(out_train[, perturbation, drop = FALSE]), perturbation),
        feat_train
      )

      test_df <- cbind(
        setNames(as.data.frame(out_test[, perturbation, drop = FALSE]), perturbation),
        feat_test
      )

      # ---- Train Model ----
      fit <- make_xgb_model(
        perturbation = perturbation,
        indx = i,
        total = total,
        dataset = train_df,
        response_cutoff = response_cutoff,
        decreasing = decreasing,
        weight_cap = weight_cap,
        nfolds = nfolds,
        nrepeats = nrepeats,
        nrounds = nrounds,
        max_depth = max_depth,
        f_subsample = f_subsample,
        min_score = min_score,
        skip_eval = skip_eval,
        shuffle = shuffle,
        n_threads = n_threads,
        xgb_params = xgb_params,
        cor_data = cor_data,
        cor_n_features = cor_n_features,
        use_gpu = use_gpu,
        gpu_id = gpu_id
      )

      # ---- Write Model Performance ----
      # metrics.json
      metrics <- list(
        jobId = job_id,
        modelKey = mk,
        perturbation = perturbation,
        mean_r = if (!is.null(fit$scores)) mean(fit$scores, na.rm = TRUE) else NA_real_,
        mean_r2 = if (!is.null(fit$scores)) mean(fit$scores^2, na.rm = TRUE) else NA_real_,
        mean_rmse = if (!is.null(fit$scores_rmse)) mean(fit$scores_rmse, na.rm = TRUE) else NA_real_,
        n_scores = if (!is.null(fit$scores)) length(fit$scores) else 0L,
        skipped = is.null(fit$model) # fit will only return a model if it passes min_score
      )
      powerup_write_json(file.path(model_out_dir, "metrics.json"), metrics)

      # ---- Make Predictions ----
      # pred_test.csv: predictions on held-out test set using the final trained model (if present).
      # pred_user.csv: predictions on user's samples using the final trained model (if present).
      # shap_user.csv: SHAP explanation values for user's sample predictions.
      if (!is.null(fit$model)) {

        # 1) Test predictions
        fit_test <- make_new_data_predictions(
          model = fit,
          name = perturbation,
          indx = i,
          total = total,
          new_data = test_df
        )

        # This added $new_data with preds/errors/shap computed from model$model + model$error_model
        pred_test_tbl <- tibble::tibble(
          cell_line = names(fit_test$new_data$predictions),
          pred = as.numeric(fit_test$new_data$predictions),
          pred_error = as.numeric(fit_test$new_data$predictions_error)
        )

        readr::write_csv(pred_test_tbl, file.path(model_out_dir, "pred_test.csv"))

        # TODO: If we want SHAP for test too, write it out:
        # shap_test_df <- fit_test$new_data$shap_values %>% tibble::rownames_to_column("cell_line")
        # readr::write_csv(shap_test_df, file.path(model_out_dir, "shap_test.csv"))

        # 2) User predictions
        fit_user <- make_new_data_predictions(
          model = fit,
          name = perturbation,
          indx = i,
          total = total,
          new_data = user_df
        )

        # This added $new_data with preds/errors/shap computed from model$model + model$error_model
        pred_user_tbl <- tibble::tibble(
          cell_line = names(fit_user$new_data$predictions),
          pred = as.numeric(fit_user$new_data$predictions),
          pred_error = as.numeric(fit_user$new_data$predictions_error)
        )

        readr::write_csv(pred_user_tbl, file.path(model_out_dir, "pred_user.csv"))

        shap_user_df <- fit_user$new_data$shap_values %>%
          as.data.frame() %>%
          tibble::rownames_to_column("cell_line")

        readr::write_csv(shap_user_df, file.path(model_out_dir, "shap_user.csv"))

        # TODO: If we also want a compact feature importance table:
        # readr::write_csv(fit_user$new_data$feature_contribution, file.path(model_out_dir, "shap_importance.csv"))

        # clean intermediates
        rm(fit_test, fit_user)
        gc()

      } else {
        # Model skipped — we still write empty artifacts
        readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_test.csv"))
        readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_user.csv"))
        readr::write_csv(tibble::tibble(), file.path(model_out_dir, "shap_user.csv"))
      }

      message(glue("[powerup][jobId={job_id}] TRAIN modelKey={mk} perturbation={perturbation} ({i} of {total}) OK skipped={is.null(fit$model)}"))

      # Clean memory between models
      rm(fit)
      gc()

    }, error = function(e) {

      # Ensure we always emit *some* artifacts for this model.
      # Keep logs very verbose so we can pinpoint the source.
      err_txt <- paste0(conditionMessage(e))
      message(glue("[powerup][jobId={job_id}] TRAIN modelKey={mk} perturbation={perturbation} ({i} of {total}) FAILED: {err_txt}"))

      # Capture a compact traceback as text (best-effort, no extra deps)
      tb <- tryCatch({
        paste(utils::capture.output(traceback(2)), collapse = "\n")
      }, error = function(.e2) {
        NA_character_
      })

      # metrics.json with error details
      metrics <- list(
        jobId = job_id,
        modelKey = mk,
        perturbation = perturbation,
        mean_r = NA_real_,
        mean_r2 = NA_real_,
        mean_rmse = NA_real_,
        n_scores = 0L,
        skipped = TRUE,
        error = list(
          message = err_txt,
          class = class(e)[1],
          traceback = tb
        )
      )
      powerup_write_json(file.path(model_out_dir, "metrics.json"), metrics)

      # Write empty artifacts matching the normal contract
      readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_test.csv"))
      readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_user.csv"))
      readr::write_csv(tibble::tibble(), file.path(model_out_dir, "shap_user.csv"))

      # Defensive cleanup
      gc()
    })
  }



  invisible(TRUE)
}

# ---- CONTRACT 3: finalize ----
#' @export
powerup_finalize <- function(shard_manifests_dir, out_aggregates_dir, job_id) {
  powerup_dir_create(out_aggregates_dir)

  # TODO: Aggregation logic:
  # - read shard manifests
  # - compile model_status.csv and a final manifest.json
  stop("powerup_finalize(): implement aggregation of shard outputs.")
}