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


.pu_write_parquet_required <- function(df, path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(glue::glue("[powerup] arrow R package is required to write parquet: {path}"))
  }
  arrow::write_parquet(df, path, compression = "zstd", use_dictionary = TRUE)
  TRUE
}


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

.pu_read_selected_samples <- function(path) {
  if (is.null(path)) return(NULL)
  if (!is.character(path) || length(path) != 1) return(NULL)

  path <- trimws(path)
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) return(NULL)

  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  ids <- unlist(strsplit(txt, "[\r\n,;]+", perl = TRUE), use.names = FALSE)
  ids <- unique(trimws(as.character(ids)))
  ids <- ids[nzchar(ids)]

  if (length(ids) < 1) return(NULL)
  ids
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


.pu_read_selected_perturbations <- function(path) {
  if (is.null(path)) return(NULL)
  if (!is.character(path) || length(path) != 1) return(NULL)

  path <- trimws(path)
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) return(NULL)

  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  ids <- unlist(strsplit(txt, "[\r\n,;]+", perl = TRUE), use.names = FALSE)
  ids <- trimws(as.character(ids))
  ids <- ids[nzchar(ids)]

  if (length(ids) < 1) return(NULL)
  ids
}

.pu_canonicalize_perturbation_id <- function(x, strip_ko_prefix = FALSE) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  if (strip_ko_prefix) {
    x <- sub("^ko_", "", x, ignore.case = TRUE)
  }
  x
}

.pu_resolve_targets_against_perturbations <- function(
  targets,
  perturbations,
  strip_ko_prefix = FALSE,
  job_id = NA_character_
) {
  if (is.null(targets) || length(targets) < 1) {
    stop(glue("[powerup][jobId={job_id}] targets resolved to empty list"))
  }

  targets <- trimws(as.character(targets))
  targets <- targets[nzchar(targets)]
  if (length(targets) < 1) {
    stop(glue("[powerup][jobId={job_id}] targets resolved to empty list after trimming"))
  }

  perturbations <- as.character(perturbations)

  available_tbl <- tibble::tibble(
    perturbation = perturbations,
    canonical = vapply(
      perturbations,
      .pu_canonicalize_perturbation_id,
      character(1),
      strip_ko_prefix = strip_ko_prefix
    )
  )

  dup_available <- available_tbl %>%
    dplyr::count(.data$canonical, name = "n") %>%
    dplyr::filter(.data$n > 1, nzchar(.data$canonical))

  if (nrow(dup_available) > 0) {
    stop(glue(
      "[powerup][jobId={job_id}] perturbation canonicalization is ambiguous for dataset columns, e.g. ",
      "{paste(head(dup_available$canonical, 10), collapse=', ')}"
    ))
  }

  wanted_tbl <- tibble::tibble(
    requested = targets,
    requestedOrder = seq_along(targets),
    canonical = vapply(
      targets,
      .pu_canonicalize_perturbation_id,
      character(1),
      strip_ko_prefix = strip_ko_prefix
    )
  ) %>%
    dplyr::group_by(.data$canonical) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  resolved_tbl <- wanted_tbl %>%
    dplyr::left_join(available_tbl, by = "canonical") %>%
    dplyr::arrange(.data$requestedOrder)

  unknown <- resolved_tbl %>%
    dplyr::filter(is.na(.data$perturbation)) %>%
    dplyr::pull(.data$requested)

  if (length(unknown) > 0) {
    message(glue(
      "[powerup][jobId={job_id}] ignoring unknown requested perturbations after exact canonical matching: ",
      "{paste(head(unknown, 25), collapse=', ')}",
      ifelse(length(unknown) > 25, " ...", "")
    ))
  }

  resolved <- resolved_tbl %>%
    dplyr::filter(!is.na(.data$perturbation)) %>%
    dplyr::pull(.data$perturbation)

  resolved <- unique(resolved)

  if (length(resolved) < 1) {
    stop(glue(
      "[powerup][jobId={job_id}] no requested perturbations matched dataset columns after exact canonical matching"
    ))
  }

  resolved
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
  selected_samples_path = NULL,
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
  if (length(common_genes) < 2) {
    stop(glue("[powerup][jobId={job_id}] Too few overlapping genes between user matrix and gene_expression: overlap={length(common_genes)}"))
  }

  perturbations <- setdiff(colnames(response_df), "cell_line")
  if (length(perturbations) < 1) stop(glue("[powerup][jobId={job_id}] response_df has no perturbation columns"))

  # ---- Optional perturbation selection (targets)
  # targets can be:
  #   - NULL: use all perturbations
  #   - single numeric (or numeric-like string): take first N perturbations (after sort)
  #   - character vector: explicit allowlist of perturbation names
  if (!is.null(targets)) {

    # Case A: numeric (or numeric-like string) => take first N
    if (is.numeric(targets) && length(targets) == 1) {
      n <- as.integer(targets)
      if (is.na(n) || n < 1) stop(glue("[powerup][jobId={job_id}] targets numeric must be >= 1"))
      perturbations <- sort(perturbations)
      if (n < length(perturbations)) {
        perturbations <- perturbations[seq_len(n)]
      }
      message(glue("[powerup][jobId={job_id}] targets numeric slice: n={n} => perturbations_used={length(perturbations)}"))
    } else if (is.character(targets) && length(targets) == 1 && grepl("^[0-9]+$", trimws(targets))) {
      n <- as.integer(trimws(targets))
      if (is.na(n) || n < 1) stop(glue("[powerup][jobId={job_id}] targets numeric-string must be >= 1"))
      perturbations <- sort(perturbations)
      if (n < length(perturbations)) {
        perturbations <- perturbations[seq_len(n)]
      }
      message(glue("[powerup][jobId={job_id}] targets numeric-string slice: n={n} => perturbations_used={length(perturbations)}"))

    # Case B: explicit allowlist vector
    } else {
      if (!is.character(targets) || length(targets) < 1) {
        stop(glue("[powerup][jobId={job_id}] targets must be NULL, a single number, or a character vector of perturbation names"))
      }

      targets_resolved <- .pu_resolve_targets_against_perturbations(
        targets = targets,
        perturbations = perturbations,
        strip_ko_prefix = identical(tolower(response_set), "crispr"),
        job_id = job_id
      )

      # Keep deterministic ordering for downstream model keys / shards
      perturbations <- sort(targets_resolved)

      requested_targets_clean <- unique(trimws(as.character(targets)))
      requested_targets_clean <- requested_targets_clean[nzchar(requested_targets_clean)]
      n_requested <- length(requested_targets_clean)
      n_matched <- length(targets_resolved)
      n_ignored <- n_requested - n_matched

      message(glue(
        "[powerup][jobId={job_id}] targets explicit allowlist resolved: ",
        "requested={n_requested} matched={n_matched} ignored={n_ignored} ",
        "stripKoPrefix={identical(tolower(response_set), 'crispr')}"
      ))

    }
  }

  # Universe = cell_line overlap between gene_expression and response
  all_ids <- intersect(gene_expression$cell_line, response_df$cell_line)
  all_ids <- sort(unique(all_ids))
  if (length(all_ids) < 3) stop(glue("[powerup][jobId={job_id}] Too few overlapping cell_line IDs between gene_expression and response: n={length(all_ids)}"))

  gene_expression_u <- gene_expression %>% dplyr::filter(.data$cell_line %in% all_ids)
  response_u        <- response_df %>%
    dplyr::filter(.data$cell_line %in% all_ids) %>%
    dplyr::select(.data$cell_line, dplyr::all_of(perturbations))

  # ---- Optional sample subset used ONLY for user-side feature selection ----
  selected_sample_ids_raw <- .pu_read_selected_samples(selected_samples_path)
  selected_sample_ids_used <- character(0)

  user_matrix_feature_source <- user_matrix

  if (!is.null(selected_sample_ids_raw) && length(selected_sample_ids_raw) > 0) {
    selected_sample_ids_used <- intersect(selected_sample_ids_raw, user_matrix$cell_line)
    selected_sample_ids_missing <- setdiff(selected_sample_ids_raw, user_matrix$cell_line)

    message(glue(
      "[powerup][jobId={job_id}] selected samples for feature selection provided requested={length(selected_sample_ids_raw)} ",
      "matched={length(selected_sample_ids_used)} missing={length(selected_sample_ids_missing)}"
    ))

    if (length(selected_sample_ids_missing) > 0) {
      message(glue(
        "[powerup][jobId={job_id}] selected samples missing from matrix: ",
        "{paste(head(selected_sample_ids_missing, 25), collapse=', ')}",
        ifelse(length(selected_sample_ids_missing) > 25, " ...", "")
      ))
    }

    if (length(selected_sample_ids_used) < 2) {
      message(glue(
        "[powerup][jobId={job_id}] selected_samples_path matched fewer than 2 samples in matrix.csv; ",
        "matched={length(selected_sample_ids_used)}. Falling back to all user samples for feature selection."
      ))
    } else {
      user_matrix_feature_source <- user_matrix %>%
        dplyr::filter(.data$cell_line %in% selected_sample_ids_used)
    }

  } else {
    message(glue("[powerup][jobId={job_id}] no selected samples file provided; using all user samples for feature selection"))
  }

  # ---- Feature selection: top variable genes on selected subset of user matrix ----
  n_features_int <- max(1L, as.integer(n_features))
  selected_genes <- .pu_select_top_var_genes(user_matrix_feature_source, common_genes, n_features_int)

  if (length(selected_genes) < 10) {
    stop(glue(
      "[powerup][jobId={job_id}] Feature selection produced too few genes: n={length(selected_genes)} ",
      "using user_samples_for_feature_selection={nrow(user_matrix_feature_source)}"
    ))
  }

  # Build feature tables
  # IMPORTANT:
  # - features_all always comes from the platform training universe
  # - features_user always keeps ALL uploaded user samples
  #   so predictions / SHAP still run on all samples
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
    selectedSamplesForFeatureSelection = list(
      provided = !is.null(selected_sample_ids_raw) && length(selected_sample_ids_raw) > 0,
      path = if (!is.null(selected_samples_path) && nzchar(selected_samples_path)) basename(selected_samples_path) else NULL,
      requested = if (!is.null(selected_sample_ids_raw)) length(selected_sample_ids_raw) else 0L,
      used = length(selected_sample_ids_used)
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


  # ------------------------------------------------------------
  # TRACE CAPTURE (logging-only)
  # Captures:
  #  - base::traceback()
  #  - sys.calls() call stack
  #  - rlang::last_trace() when available
  # ------------------------------------------------------------
  .pu_capture_traces <- function() {
    out <- list(
      traceback = NA_character_,
      calls = NA_character_,
      rlang_last_trace = NA_character_
    )

    out$traceback <- tryCatch(
      paste(utils::capture.output(base::traceback()), collapse = "\n"),
      error = function(e) NA_character_
    )

    out$calls <- tryCatch({
      calls <- sys.calls()
      paste(vapply(calls, function(x) paste(deparse(x), collapse = ""), character(1)), collapse = "\n")
    }, error = function(e) NA_character_)

    out$rlang_last_trace <- tryCatch({
      if (!requireNamespace("rlang", quietly = TRUE)) return(NA_character_)
      lt <- rlang::last_trace()
      paste(utils::capture.output(print(lt)), collapse = "\n")
    }, error = function(e) NA_character_)

    out
  }

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (is.character(x) && !nzchar(x))) y else x
  



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

  # Sanity check: all model_keys must be present in perturbations.csv mapping
  if (nrow(key_map) != length(model_keys)) {
    missing <- setdiff(model_keys, key_map$modelKey)
    stop(glue("Some modelKeys missing in perturbations.csv: {paste(missing, collapse=', ')}"))
  }

  # Reorder key_map to match order of model_keys input (for deterministic processing downstream)
  key_map <- key_map %>%
    mutate(.order = match(.data$modelKey, model_keys)) %>%
    arrange(.data$.order) %>%
    select(-.data$.order)

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
    mk <- key_map$modelKey[[i]]
    perturbation <- key_map$perturbation[[i]]

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
        mean_r2 = if (!is.null(fit$scores_R2)) mean(fit$scores_R2, na.rm = TRUE) else NA_real_,
        mean_rmse = if (!is.null(fit$scores_rmse)) mean(fit$scores_rmse, na.rm = TRUE) else NA_real_,
        n_scores = if (!is.null(fit$scores)) length(fit$scores) else 0L,
        skipped = is.null(fit$model)
      )
      powerup_write_json(file.path(model_out_dir, "metrics.json"), metrics)

      # ---- Make Predictions ----
      # pred_test.csv: predictions on held-out test set using the final trained model (if present).
      # pred_user.csv: predictions on user's samples using the final trained model (if present).
      # shap_user.csv: SHAP explanation values for user's sample predictions.
      if (!is.null(fit$model)) {

        # -----------------------------
        # Stage-tagged prediction calls
        # (logging-only; no behavior changes)
        # -----------------------------
        stage <- "INIT"

        fit_test <- NULL
        fit_user <- NULL

        # 1) Test predictions
        stage <- "PREDICT_TEST: make_new_data_predictions(test_df)"
        fit_test <- tryCatch({
          make_new_data_predictions(
            model = fit,
            name = perturbation,
            indx = i,
            total = total,
            new_data = test_df
          )
        }, error = function(e) {
          tr <- .pu_capture_traces()
          msg <- glue("[powerup][jobId={job_id}] {stage} FAILED modelKey={mk} perturbation={perturbation}: {conditionMessage(e)}")
          message(msg)
          if (nzchar(tr$traceback %||% "")) message(glue("[powerup][jobId={job_id}] {stage} traceback:\n{tr$traceback}"))
          if (nzchar(tr$calls %||% "")) message(glue("[powerup][jobId={job_id}] {stage} calls:\n{tr$calls}"))
          if (nzchar(tr$rlang_last_trace %||% "")) message(glue("[powerup][jobId={job_id}] {stage} rlang::last_trace:\n{tr$rlang_last_trace}"))
          stop(e)
        })

        pred_test_tbl <- tibble::tibble(
          cell_line = names(fit_test$new_data$pred_mean),
          pred = as.numeric(fit_test$new_data$pred_mean),                 # legacy-compatible alias
          pred_error = as.numeric(fit_test$new_data$pred_sd),             # legacy-compatible alias
          pred_mean = as.numeric(fit_test$new_data$pred_mean),
          pred_var = as.numeric(fit_test$new_data$pred_var),
          pred_sd = as.numeric(fit_test$new_data$pred_sd),
          pred_log_var = as.numeric(fit_test$new_data$pred_log_var),
          pred_pi_lower_95 = as.numeric(fit_test$new_data$pred_pi_lower_95),
          pred_pi_upper_95 = as.numeric(fit_test$new_data$pred_pi_upper_95),
          response_cutoff = as.numeric(fit_test$new_data$response_cutoff),
          decreasing = as.logical(fit_test$new_data$decreasing),
          prob_below_cutoff = as.numeric(fit_test$new_data$prob_below_cutoff),
          prob_above_cutoff = as.numeric(fit_test$new_data$prob_above_cutoff),
          prob_target_event = as.numeric(fit_test$new_data$prob_target_event),
          target_event_definition = as.character(fit_test$new_data$target_event_definition)
        )
        readr::write_csv(pred_test_tbl, file.path(model_out_dir, "pred_test.csv"))

        # 2) User predictions
        stage <- "PREDICT_USER: make_new_data_predictions(user_df)"
        fit_user <- tryCatch({
          make_new_data_predictions(
            model = fit,
            name = perturbation,
            indx = i,
            total = total,
            new_data = user_df
          )
        }, error = function(e) {
          tr <- .pu_capture_traces()
          msg <- glue("[powerup][jobId={job_id}] {stage} FAILED modelKey={mk} perturbation={perturbation}: {conditionMessage(e)}")
          message(msg)
          if (nzchar(tr$traceback %||% "")) message(glue("[powerup][jobId={job_id}] {stage} traceback:\n{tr$traceback}"))
          if (nzchar(tr$calls %||% "")) message(glue("[powerup][jobId={job_id}] {stage} calls:\n{tr$calls}"))
          if (nzchar(tr$rlang_last_trace %||% "")) message(glue("[powerup][jobId={job_id}] {stage} rlang::last_trace:\n{tr$rlang_last_trace}"))
          stop(e)
        })

        pred_user_tbl <- tibble::tibble(
          cell_line = names(fit_user$new_data$pred_mean),
          pred = as.numeric(fit_user$new_data$pred_mean),                 # legacy-compatible alias
          pred_error = as.numeric(fit_user$new_data$pred_sd),             # legacy-compatible alias
          pred_mean = as.numeric(fit_user$new_data$pred_mean),
          pred_var = as.numeric(fit_user$new_data$pred_var),
          pred_sd = as.numeric(fit_user$new_data$pred_sd),
          pred_log_var = as.numeric(fit_user$new_data$pred_log_var),
          pred_pi_lower_95 = as.numeric(fit_user$new_data$pred_pi_lower_95),
          pred_pi_upper_95 = as.numeric(fit_user$new_data$pred_pi_upper_95),
          response_cutoff = as.numeric(fit_user$new_data$response_cutoff),
          decreasing = as.logical(fit_user$new_data$decreasing),
          prob_below_cutoff = as.numeric(fit_user$new_data$prob_below_cutoff),
          prob_above_cutoff = as.numeric(fit_user$new_data$prob_above_cutoff),
          prob_target_event = as.numeric(fit_user$new_data$prob_target_event),
          target_event_definition = as.character(fit_user$new_data$target_event_definition)
        )
        readr::write_csv(pred_user_tbl, file.path(model_out_dir, "pred_user.csv"))


        # 3) SHAP extraction
        stage <- "SHAP_USER: extract fit_user$new_data$shap_values"
        shap_user_df <- tryCatch({
          fit_user$new_data$shap_values %>%
            as.data.frame() %>%
            tibble::rownames_to_column("cell_line")
        }, error = function(e) {
          tr <- .pu_capture_traces()
          msg <- glue("[powerup][jobId={job_id}] {stage} FAILED modelKey={mk} perturbation={perturbation}: {conditionMessage(e)}")
          message(msg)
          if (nzchar(tr$traceback %||% "")) message(glue("[powerup][jobId={job_id}] {stage} traceback:\n{tr$traceback}"))
          if (nzchar(tr$calls %||% "")) message(glue("[powerup][jobId={job_id}] {stage} calls:\n{tr$calls}"))
          if (nzchar(tr$rlang_last_trace %||% "")) message(glue("[powerup][jobId={job_id}] {stage} rlang::last_trace:\n{tr$rlang_last_trace}"))
          stop(e)
        })

        #readr::write_csv(shap_user_df, file.path(model_out_dir, "shap_user.csv"))
        .pu_write_parquet_required(shap_user_df, file.path(model_out_dir, "shap_user.parquet"))

        rm(fit_test, fit_user)
        rm(shap_user_df)
        gc()


        # 4) TRAINING SHAP extraction (same contract as SHAP_USER)
        stage <- "SHAP_TRAIN: extract fit$shap_values"
        shap_train_df <- tryCatch({
          fit$shap_values %>%
            as.data.frame() %>%
            tibble::rownames_to_column("cell_line")
        }, error = function(e) {
          tr <- .pu_capture_traces()
          msg <- glue("[powerup][jobId={job_id}] {stage} FAILED modelKey={mk} perturbation={perturbation}: {conditionMessage(e)}")
          message(msg)
          if (nzchar(tr$traceback %||% "")) message(glue("[powerup][jobId={job_id}] {stage} traceback:\n{tr$traceback}"))
          if (nzchar(tr$calls %||% "")) message(glue("[powerup][jobId={job_id}] {stage} calls:\n{tr$calls}"))
          if (nzchar(tr$rlang_last_trace %||% "")) message(glue("[powerup][jobId={job_id}] {stage} rlang::last_trace:\n{tr$rlang_last_trace}"))
          stop(e)
        })

        # readr::write_csv(shap_train_df, file.path(model_out_dir, "shap_train.csv"))
        .pu_write_parquet_required(shap_train_df, file.path(model_out_dir, "shap_train.parquet"))

        rm(shap_train_df)
        gc()

      } else {
        readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_test.csv"))
        readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_user.csv"))
        .pu_write_parquet_required(tibble::tibble(), file.path(model_out_dir, "shap_user.parquet"))
        .pu_write_parquet_required(tibble::tibble(), file.path(model_out_dir, "shap_train.parquet"))
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
      tr <- .pu_capture_traces()
      tb <- tr$traceback
      calls_txt <- tr$calls
      rlang_trace_txt <- tr$rlang_last_trace

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
          traceback = tb,
          calls = calls_txt,
          rlang_last_trace = rlang_trace_txt
        )
      )
      powerup_write_json(file.path(model_out_dir, "metrics.json"), metrics)

      # Write empty artifacts matching the normal contract
      readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_test.csv"))
      readr::write_csv(tibble::tibble(), file.path(model_out_dir, "pred_user.csv"))
      .pu_write_parquet_required(tibble::tibble(), file.path(model_out_dir, "shap_user.parquet"))
      .pu_write_parquet_required(tibble::tibble(), file.path(model_out_dir, "shap_train.parquet"))

      # Defensive cleanup
      gc()
    })
  }



  invisible(TRUE)
}

# ---- CONTRACT 3: finalize ----
# Aggregates model-level metrics.json into a single model_status.csv.
# Shard-agnostic and local-only.
# Expects:
#   - out_preprocess_dir contains perturbations.csv (written by powerup_preprocess)
#   - models_dir contains models/<modelKey>/metrics.json (may be missing for some models)
#' @export
powerup_finalize <- function(out_preprocess_dir, models_dir, out_aggregates_dir, job_id) {
  powerup_dir_create(out_aggregates_dir)

  .pu_assert_file_exists(out_preprocess_dir, "out_preprocess_dir")
  if (!dir.exists(out_preprocess_dir)) stop(glue("out_preprocess_dir is not a directory: {out_preprocess_dir}"))

  .pu_assert_file_exists(models_dir, "models_dir")
  if (!dir.exists(models_dir)) stop(glue("models_dir is not a directory: {models_dir}"))

  perturbations_path <- file.path(out_preprocess_dir, "perturbations.csv")
  .pu_assert_file_exists(perturbations_path, "perturbations.csv (in out_preprocess_dir)")

  message(glue(
    "[powerup][jobId={job_id}] FINALIZE start preprocess_dir={out_preprocess_dir} models_dir={models_dir} out={out_aggregates_dir}"
  ))

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

  # ---- Source of truth: expected models + perturbation names ----
  pert_tbl <- readr::read_csv(perturbations_path, show_col_types = FALSE, progress = FALSE)
  if (!("modelKey" %in% names(pert_tbl))) stop("perturbations.csv must contain column: modelKey")
  if (!("perturbation" %in% names(pert_tbl))) stop("perturbations.csv must contain column: perturbation")

  pert_tbl <- pert_tbl %>%
    dplyr::mutate(
      modelKey = as.character(.data$modelKey),
      perturbation = as.character(.data$perturbation)
    ) %>%
    dplyr::filter(!is.na(.data$modelKey), trimws(.data$modelKey) != "") %>%
    dplyr::distinct(.data$modelKey, .keep_all = TRUE) %>%
    dplyr::arrange(.data$modelKey)

  if (nrow(pert_tbl) < 1) stop(glue("[powerup][jobId={job_id}] perturbations.csv has 0 valid modelKey rows"))

  message(glue("[powerup][jobId={job_id}] FINALIZE expected models={nrow(pert_tbl)}"))

  # Deterministic platform paths (still OK in fully local runs as "intended locations")
  gcs_prefix <- glue("jobs/{job_id}/output/models")

  # ---- Discover local metrics.json ----
  metric_files <- list.files(models_dir, pattern = "^metrics\\.json$", recursive = TRUE, full.names = TRUE)
  metric_files <- sort(unique(metric_files))

  metrics_map <- tibble::tibble(
    metricsFile = metric_files,
    modelKey = if (length(metric_files) > 0) basename(dirname(metric_files)) else character(0)
  ) %>%
    dplyr::filter(!is.na(.data$modelKey), trimws(.data$modelKey) != "") %>%
    dplyr::arrange(.data$modelKey, .data$metricsFile) %>%
    dplyr::group_by(.data$modelKey) %>%
    dplyr::summarise(metricsFile = dplyr::first(.data$metricsFile), .groups = "drop")

  message(glue("[powerup][jobId={job_id}] FINALIZE discovered metrics.json={nrow(metrics_map)}"))

  parse_metrics <- function(path, fallback_model_key, fallback_perturbation) {
    jobId_m <- job_id
    perturbation <- fallback_perturbation
    mean_r <- NA_real_
    mean_r2 <- NA_real_
    mean_rmse <- NA_real_
    n_scores <- 0L
    skipped <- TRUE
    errorType <- NA_character_
    errorMessage <- NA_character_

    txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
    if (is.null(txt) || length(txt) == 0) {
      errorType <- "LocalReadError"
      errorMessage <- glue("Could not read file: {path}")
    } else {
      obj <- tryCatch(
        jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = TRUE),
        error = function(e) {
          errorType <<- "JsonParseError"
          errorMessage <<- conditionMessage(e)
          NULL
        }
      )

      if (!is.null(obj)) {
        jobId_m <- as.character(obj$jobId %||% job_id)
        perturbation <- as.character(obj$perturbation %||% fallback_perturbation)

        mean_r <- suppressWarnings(as.numeric(obj$mean_r %||% NA_real_))
        mean_r2 <- suppressWarnings(as.numeric(obj$mean_r2 %||% NA_real_))
        mean_rmse <- suppressWarnings(as.numeric(obj$mean_rmse %||% NA_real_))
        n_scores <- suppressWarnings(as.integer(obj$n_scores %||% 0L))
        skipped <- isTRUE(obj$skipped %||% FALSE)

        if (!is.null(obj$error)) {
          errorType <- as.character(obj$error$class %||% NA_character_)
          errorMessage <- as.character(obj$error$message %||% NA_character_)
        }
      }
    }

    status <- "OK"
    if (!is.na(errorType) || !is.na(errorMessage)) {
      status <- "FAILED"
    } else if (isTRUE(skipped)) {
      status <- "SKIPPED"
    }

    model_dir_gcs <- glue("{gcs_prefix}/{fallback_model_key}")

    tibble::tibble(
      jobId = jobId_m,
      modelKey = fallback_model_key,
      perturbation = perturbation,
      status = status,
      skipped = isTRUE(skipped),
      mean_r = mean_r,
      mean_r2 = mean_r2,
      mean_rmse = mean_rmse,
      n_scores = n_scores,
      errorType = errorType,
      errorMessage = errorMessage,
      metricsPath = glue("{model_dir_gcs}/metrics.json"),
      predTestPath = glue("{model_dir_gcs}/pred_test.csv"),
      predUserPath = glue("{model_dir_gcs}/pred_user.csv"),
      shapUserPath = glue("{model_dir_gcs}/shap_user.parquet"),
      shapTrainPath = glue("{model_dir_gcs}/shap_train.parquet")
    )
  }

  base <- pert_tbl %>%
    dplyr::select(.data$modelKey, .data$perturbation) %>%
    dplyr::left_join(metrics_map, by = "modelKey")

  rows <- vector("list", nrow(base))

  for (i in seq_len(nrow(base))) {
    mk <- base$modelKey[[i]]
    pert <- base$perturbation[[i]]
    mf <- base$metricsFile[[i]]

    model_dir_gcs <- glue("{gcs_prefix}/{mk}")

    if (is.na(mf) || !nzchar(mf) || !file.exists(mf)) {
      rows[[i]] <- tibble::tibble(
        jobId = job_id,
        modelKey = mk,
        perturbation = pert,
        status = "MISSING_METRICS",
        skipped = TRUE,
        mean_r = NA_real_,
        mean_r2 = NA_real_,
        mean_rmse = NA_real_,
        n_scores = 0L,
        errorType = "MissingMetrics",
        errorMessage = "metrics.json not found locally for this modelKey",
        metricsPath = glue("{model_dir_gcs}/metrics.json"),
        predTestPath = glue("{model_dir_gcs}/pred_test.csv"),
        predUserPath = glue("{model_dir_gcs}/pred_user.csv"),
        shapUserPath = glue("{model_dir_gcs}/shap_user.parquet"),
        shapTrainPath = glue("{model_dir_gcs}/shap_train.parquet")
      )
    } else {
      rows[[i]] <- parse_metrics(mf, fallback_model_key = mk, fallback_perturbation = pert)
    }
  }

  model_status <- dplyr::bind_rows(rows) %>%
    dplyr::arrange(.data$modelKey)

  out_csv <- file.path(out_aggregates_dir, "model_status.csv")
  readr::write_csv(model_status, out_csv)

  total <- nrow(model_status)
  n_ok <- sum(model_status$status == "OK")
  n_skipped <- sum(model_status$status == "SKIPPED")
  n_failed <- sum(model_status$status == "FAILED")
  n_missing <- sum(model_status$status == "MISSING_METRICS")

  powerup_write_json(file.path(out_aggregates_dir, "manifest.json"), list(
    schemaVersion = 1,
    jobId = job_id,
    createdAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
    counts = list(
      totalModels = total,
      ok = n_ok,
      skipped = n_skipped,
      failed = n_failed,
      missingMetrics = n_missing
    ),
    artifacts = list(
      modelStatusCsv = "model_status.csv"
    )
  ))

  message(glue("[powerup][jobId={job_id}] FINALIZE done total={total} ok={n_ok} skipped={n_skipped} failed={n_failed} missingMetrics={n_missing}"))
  invisible(TRUE)
}




#' Aggregate user predictions across passing models (local runs)
#'
#' This is a convenience helper for users running PowerUp locally. It mirrors the
#' platform FINALIZE aggregation behavior by reading model_status.csv and
#' combining per-model pred_user.csv into a single long table, plus a small JSON index.
#'
#' Passing models are defined as:
#'   - status == "OK"
#'   - skipped == FALSE
#'
#' The function attempts to resolve pred_user.csv paths robustly:
#'   1) If predUserPath exists locally as-is, it is used.
#'   2) If predUserPath looks like a platform path (e.g., jobs/<jobId>/output/models/<modelKey>/pred_user.csv),
#'      it is mapped to file.path(models_dir, <modelKey>, "pred_user.csv").
#'   3) Otherwise it falls back to file.path(models_dir, <modelKey>, "pred_user.csv").
#'
#' Outputs (written to out_aggregates_dir):
#'   - user_predictions_long.csv
#'   - user_predictions_index.json
#'
#' Determinism:
#'   - Models processed in sorted(modelKey) order.
#'   - Rows sorted by cell_line within each model before writing.
#'
#' @param out_aggregates_dir Directory containing model_status.csv and where aggregated output is written.
#' @param models_dir Directory containing per-model outputs under models_dir/<modelKey>/pred_user.csv.
#' @param model_status_path Optional explicit path to model_status.csv.
#'        Defaults to file.path(out_aggregates_dir, "model_status.csv").
#' @param output_filename Name of the aggregated CSV to write in out_aggregates_dir.
#' @param index_filename Name of the index JSON to write in out_aggregates_dir.
#' @param verbose If TRUE, prints progress + warnings.
#'
#' @return A list with counts and output paths. (Invisible)
#' @export
powerup_aggregate_user_predictions <- function(
  out_aggregates_dir,
  models_dir,
  model_status_path = file.path(out_aggregates_dir, "model_status.csv"),
  output_filename = "user_predictions_long.csv",
  index_filename = "user_predictions_index.json",
  verbose = TRUE
) {
  powerup_dir_create(out_aggregates_dir)

  .pu_assert_file_exists(out_aggregates_dir, "out_aggregates_dir")
  if (!dir.exists(out_aggregates_dir)) stop(glue("out_aggregates_dir is not a directory: {out_aggregates_dir}"))

  .pu_assert_file_exists(models_dir, "models_dir")
  if (!dir.exists(models_dir)) stop(glue("models_dir is not a directory: {models_dir}"))

  .pu_assert_file_exists(model_status_path, "model_status_path")

  if (!is.character(output_filename) || length(output_filename) != 1 || !nzchar(output_filename)) {
    stop("output_filename must be a non-empty string")
  }
  if (!is.character(index_filename) || length(index_filename) != 1 || !nzchar(index_filename)) {
    stop("index_filename must be a non-empty string")
  }

  # ---- Read model_status.csv ----
  ms <- readr::read_csv(model_status_path, show_col_types = FALSE, progress = FALSE)

  required_cols <- c("jobId", "modelKey", "perturbation", "status", "skipped", "predUserPath")
  missing <- setdiff(required_cols, names(ms))
  if (length(missing) > 0) {
    stop(glue(
      "model_status.csv is missing required columns: {paste(missing, collapse = ', ')}"
    ))
  }

  # Normalize types
  ms <- ms %>%
    dplyr::mutate(
      jobId = as.character(.data$jobId),
      modelKey = as.character(.data$modelKey),
      perturbation = as.character(.data$perturbation),
      status = as.character(.data$status),
      skipped = as.logical(.data$skipped),
      predUserPath = as.character(.data$predUserPath),
      mean_r = if ("mean_r" %in% names(ms)) suppressWarnings(as.numeric(.data$mean_r)) else NA_real_,
      mean_rmse = if ("mean_rmse" %in% names(ms)) suppressWarnings(as.numeric(.data$mean_rmse)) else NA_real_
    ) %>%
    dplyr::filter(!is.na(.data$modelKey), trimws(.data$modelKey) != "") %>%
    dplyr::arrange(.data$modelKey)

  passed <- ms %>%
    dplyr::filter(.data$status == "OK", isFALSE(.data$skipped)) %>%
    dplyr::arrange(.data$modelKey)

  out_csv_path <- file.path(out_aggregates_dir, output_filename)
  out_index_path <- file.path(out_aggregates_dir, index_filename)

  # ---- Path resolver ----
  resolve_pred_user_path <- function(predUserPath, modelKey) {
    p <- predUserPath
    if (!is.na(p) && nzchar(p) && file.exists(p)) return(p)

    # platform-like relative path:
    # jobs/<jobId>/output/models/<modelKey>/pred_user.csv
    if (!is.na(p) && nzchar(p) && grepl("^jobs/.+/output/models/.+/pred_user\\.csv$", p)) {
      return(file.path(models_dir, modelKey, "pred_user.csv"))
    }

    # fallback
    file.path(models_dir, modelKey, "pred_user.csv")
  }

  # ---- Prepare streaming CSV output ----
  out_con <- file(out_csv_path, open = "wt")
  on.exit(close(out_con), add = TRUE)

  header <- c(
    "jobId",
    "modelKey",
    "perturbation",
    "cell_line",
    "pred",
    "pred_error",
    "metrics_mean_r",
    "metrics_mean_rmse"
  )
  writeLines(paste(header, collapse = ","), con = out_con)

  # For index.json
  index_models <- list()

  n_models_with_preds <- 0L
  user_samples_seen <- new.env(parent = emptyenv())  # set-like
  n_written <- 0L

  # Minimal CSV escaping for strings
  csv_escape <- function(x) {
    if (is.na(x)) return("")
    s <- as.character(x)
    if (grepl('[,"\r\n]', s)) {
      s <- gsub('"', '""', s, fixed = TRUE)
      return(paste0('"', s, '"'))
    }
    s
  }

  # If no passing models, write empty outputs deterministically
  if (nrow(passed) < 1) {
    # CSV already has header; nothing else to write
    idx <- list(
      schemaVersion = 1,
      generatedAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
      counts = list(
        nModelsTotal = nrow(ms),
        nModelsPassed = 0L,
        nModelsWithUserPreds = 0L,
        nUserSamples = 0L
      ),
      artifacts = list(
        user_predictions_long = basename(out_csv_path),
        user_predictions_index = basename(out_index_path)
      ),
      models = list()
    )
    jsonlite::write_json(idx, out_index_path, auto_unbox = TRUE, pretty = TRUE)
    if (isTRUE(verbose)) {
      message(glue("[powerup] aggregate_user_predictions: no passing models; wrote empty {out_csv_path} and {out_index_path}"))
    }
    return(invisible(list(
      ok = TRUE,
      outputs = list(csv = out_csv_path, index = out_index_path),
      counts = idx$counts
    )))
  }

  # ---- Process each passing model deterministically ----
  for (i in seq_len(nrow(passed))) {
    mk <- passed$modelKey[[i]]
    jobId_i <- passed$jobId[[i]]
    pert <- passed$perturbation[[i]]
    mean_r <- passed$mean_r[[i]]
    mean_rmse <- passed$mean_rmse[[i]]
    pred_path <- resolve_pred_user_path(passed$predUserPath[[i]], mk)

    # Default index record (updated below)
    idx_rec <- list(
      modelKey = mk,
      perturbation = pert,
      passed = TRUE,
      hasUserPreds = FALSE,
      predUserPath = pred_path
    )

    if (!file.exists(pred_path)) {
      if (isTRUE(verbose)) message(glue("[powerup] aggregate_user_predictions: missing pred_user.csv; modelKey={mk} path={pred_path}"))
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    dt <- tryCatch(
      data.table::fread(pred_path, data.table = FALSE, showProgress = FALSE),
      error = function(e) e
    )

    if (inherits(dt, "error") || is.null(dt) || nrow(dt) < 1) {
      if (isTRUE(verbose)) message(glue("[powerup] aggregate_user_predictions: empty/invalid pred_user.csv; modelKey={mk} path={pred_path}"))
      if (inherits(dt, "error")) {
        idx_rec$error <- list(type = class(dt)[1], message = conditionMessage(dt))
      }
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    # Normalize ID col
    if (!("cell_line" %in% colnames(dt))) {
      colnames(dt)[1] <- "cell_line"
    }
    dt$cell_line <- as.character(dt$cell_line)

    fields <- colnames(dt)

    # Choose prediction column
    pred_col <- if ("pred" %in% fields) {
      "pred"
    } else if ("prediction" %in% fields) {
      "prediction"
    } else {
      setdiff(fields, "cell_line")[1]
    }

    if (is.na(pred_col) || !nzchar(pred_col) || !(pred_col %in% fields)) {
      if (isTRUE(verbose)) message(glue("[powerup] aggregate_user_predictions: cannot find pred column; modelKey={mk} path={pred_path}"))
      idx_rec$error <- list(type = "SchemaError", message = "Could not determine prediction column")
      idx_rec$predUserSchema <- list(
        cell_col = "cell_line",
        pred_col = NA,
        pred_error_col = NA,
        columns = fields
      )
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    # Optional pred_error column
    pred_err_col <- if ("pred_error" %in% fields) {
      "pred_error"
    } else if ("se" %in% fields) {
      "se"
    } else {
      NA_character_
    }

    # Coerce numeric safely
    suppressWarnings(dt[[pred_col]] <- as.numeric(dt[[pred_col]]))
    if (!is.na(pred_err_col) && (pred_err_col %in% fields)) {
      suppressWarnings(dt[[pred_err_col]] <- as.numeric(dt[[pred_err_col]]))
    } else {
      dt$pred_error <- NA_real_
      pred_err_col <- "pred_error"
    }

    # Build standardized output + deterministic ordering
    out_tbl <- dt %>%
      dplyr::transmute(
        jobId = jobId_i,
        modelKey = mk,
        perturbation = pert,
        cell_line = .data$cell_line,
        pred = .data[[pred_col]],
        pred_error = .data[[pred_err_col]],
        metrics_mean_r = mean_r,
        metrics_mean_rmse = mean_rmse
      ) %>%
      dplyr::filter(!is.na(.data$cell_line), trimws(.data$cell_line) != "") %>%
      dplyr::arrange(.data$cell_line)

    if (nrow(out_tbl) < 1) {
      idx_rec$predUserSchema <- list(
        cell_col = "cell_line",
        pred_col = pred_col,
        pred_error_col = ifelse(is.na(pred_err_col), NA, pred_err_col),
        columns = fields
      )
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    # Track unique user samples
    for (cl in out_tbl$cell_line) user_samples_seen[[cl]] <- TRUE

    # Stream write rows
    for (r in seq_len(nrow(out_tbl))) {
      row <- out_tbl[r, ]
      line <- paste(
        csv_escape(row$jobId),
        csv_escape(row$modelKey),
        csv_escape(row$perturbation),
        csv_escape(row$cell_line),
        ifelse(is.na(row$pred), "", format(row$pred, scientific = FALSE)),
        ifelse(is.na(row$pred_error), "", format(row$pred_error, scientific = FALSE)),
        ifelse(is.na(row$metrics_mean_r), "", format(row$metrics_mean_r, scientific = FALSE)),
        ifelse(is.na(row$metrics_mean_rmse), "", format(row$metrics_mean_rmse, scientific = FALSE)),
        sep = ","
      )
      writeLines(line, con = out_con)
      n_written <- n_written + 1L
    }

    n_models_with_preds <- n_models_with_preds + 1L

    idx_rec$hasUserPreds <- TRUE
    idx_rec$mean_r <- mean_r
    idx_rec$mean_rmse <- mean_rmse
    idx_rec$predUserSchema <- list(
      cell_col = "cell_line",
      pred_col = pred_col,
      pred_error_col = pred_err_col,
      columns = fields
    )

    index_models[[length(index_models) + 1L]] <- idx_rec

    if (isTRUE(verbose)) message(glue("[powerup] aggregate_user_predictions: aggregated modelKey={mk} rows={nrow(out_tbl)}"))
  }

  n_user_samples <- length(ls(user_samples_seen, all.names = TRUE))

  # ---- Write index JSON ----
  idx <- list(
    schemaVersion = 1,
    generatedAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
    counts = list(
      nModelsTotal = nrow(ms),
      nModelsPassed = nrow(passed),
      nModelsWithUserPreds = n_models_with_preds,
      nUserSamples = n_user_samples
    ),
    artifacts = list(
      user_predictions_long = basename(out_csv_path),
      user_predictions_index = basename(out_index_path)
    ),
    # Deterministic ordering already enforced by passed %>% arrange(modelKey)
    models = index_models
  )
  jsonlite::write_json(idx, out_index_path, auto_unbox = TRUE, pretty = TRUE)

  if (isTRUE(verbose)) {
    message(glue(
      "[powerup] aggregate_user_predictions: wrote {out_csv_path} rows={n_written} models_passed={nrow(passed)} models_with_user_preds={n_models_with_preds} user_samples={n_user_samples}"
    ))
    message(glue("[powerup] aggregate_user_predictions: wrote {out_index_path} models_indexed={length(index_models)}"))
  }

  invisible(list(
    ok = TRUE,
    outputs = list(csv = out_csv_path, index = out_index_path),
    counts = idx$counts
  ))
}




#' Aggregate TRAINING SHAP values across passing models (local runs)
#'
#' Reads per-model shap_train.parquet (WIDE) and writes a LONG parquet "dataset"
#' as multiple part files + a JSON index.
#'
#' Passing models:
#'   - status == "OK"
#'   - skipped == FALSE
#'
#' Determinism:
#'   - Models processed in sorted(modelKey) order.
#'   - cell_line order preserved within each parquet batch (arrow batch order).
#'   - Features ordered with "(Intercept)" first (if present), then alphabetical.
#'
#' Filtering:
#'   - Keeps ONLY abs(shap_value) > 0 (non-zero)
#'   - Excludes "(Intercept)" from output rows (consistent with worker finalize melt)
#'
#' Outputs:
#'   - A directory under out_aggregates_dir: <output_basename>/
#'       - part-00001.parquet, part-00002.parquet, ...
#'   - training_shap_index.json
#'
#' @export
powerup_aggregate_training_shap <- function(
  out_aggregates_dir,
  models_dir,
  model_status_path = file.path(out_aggregates_dir, "model_status.csv"),
  output_filename = "training_shap_long",   # will be a DIRECTORY, not a single file
  index_filename = "training_shap_index.json",
  verbose = TRUE,
  part_rows = 500000L  # rows per parquet part (tune as needed)
) {
  powerup_dir_create(out_aggregates_dir)

  .pu_assert_file_exists(out_aggregates_dir, "out_aggregates_dir")
  if (!dir.exists(out_aggregates_dir)) stop(glue("out_aggregates_dir is not a directory: {out_aggregates_dir}"))

  .pu_assert_file_exists(models_dir, "models_dir")
  if (!dir.exists(models_dir)) stop(glue("models_dir is not a directory: {models_dir}"))

  .pu_assert_file_exists(model_status_path, "model_status_path")

  if (!is.character(output_filename) || length(output_filename) != 1 || !nzchar(output_filename)) {
    stop("output_filename must be a non-empty string (directory basename)")
  }
  if (!is.character(index_filename) || length(index_filename) != 1 || !nzchar(index_filename)) {
    stop("index_filename must be a non-empty string")
  }

  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(glue("[powerup] arrow R package is required to read parquet shap_train.parquet"))
  }

  # ---- Read model_status.csv ----
  ms <- readr::read_csv(model_status_path, show_col_types = FALSE, progress = FALSE)

  required_cols <- c("jobId", "modelKey", "perturbation", "status", "skipped")
  missing <- setdiff(required_cols, names(ms))
  if (length(missing) > 0) {
    stop(glue("model_status.csv is missing required columns: {paste(missing, collapse = ', ')}"))
  }

  ms <- ms %>%
    dplyr::mutate(
      jobId = as.character(.data$jobId),
      modelKey = as.character(.data$modelKey),
      perturbation = as.character(.data$perturbation),
      status = as.character(.data$status),
      skipped = as.logical(.data$skipped),
      shapTrainPath = if ("shapTrainPath" %in% names(ms)) as.character(.data$shapTrainPath) else NA_character_
    ) %>%
    dplyr::filter(!is.na(.data$modelKey), trimws(.data$modelKey) != "") %>%
    dplyr::arrange(.data$modelKey)

  passed <- ms %>%
    dplyr::filter(.data$status == "OK", isFALSE(.data$skipped)) %>%
    dplyr::arrange(.data$modelKey)

  # ---- Output locations ----
  out_dir <- file.path(out_aggregates_dir, output_filename)
  powerup_dir_create(out_dir)
  out_index_path <- file.path(out_aggregates_dir, index_filename)

  # Clean existing part files deterministically (local helper; safe)
  old_parts <- list.files(out_dir, pattern = "^part-[0-9]+\\.parquet$", full.names = TRUE)
  if (length(old_parts) > 0) file.remove(old_parts)

  # ---- Path resolver: expects shap_train.parquet now ----
  resolve_shap_train_path <- function(shapTrainPath, modelKey) {
    p <- shapTrainPath
    if (!is.na(p) && nzchar(p) && file.exists(p)) return(p)

    # platform-like relative path:
    # jobs/<jobId>/output/models/<modelKey>/shap_train.parquet
    if (!is.na(p) && nzchar(p) && grepl("^jobs/.+/output/models/.+/shap_train\\.parquet$", p)) {
      return(file.path(models_dir, modelKey, "shap_train.parquet"))
    }

    # fallback
    file.path(models_dir, modelKey, "shap_train.parquet")
  }

  # ---- Index accumulators ----
  index_models <- list()
  n_models_with_training_shap <- 0L
  n_rows_total <- 0L
  train_samples_seen <- new.env(parent = emptyenv())

  # Buffer for writing part files
  buf <- NULL
  buf_n <- 0L
  part_i <- 0L

  # Schema for long output
  long_schema <- arrow::schema(
    jobId = arrow::utf8(),
    modelKey = arrow::utf8(),
    perturbation = arrow::utf8(),
    cell_line = arrow::utf8(),
    feature = arrow::utf8(),
    shap_value = arrow::float64()
  )

  flush_part <- function() {
    if (is.null(buf) || buf_n < 1L) return(invisible(FALSE))
    part_i <<- part_i + 1L
    part_path <- file.path(out_dir, sprintf("part-%05d.parquet", part_i))

    # write_parquet expects an Arrow Table or data.frame; we use Arrow Table
    arrow::write_parquet(buf, part_path, compression = "zstd", use_dictionary = TRUE)

    if (isTRUE(verbose)) {
      message(glue("[powerup] aggregate_training_shap: wrote {basename(part_path)} rows={buf_n}"))
    }

    n_rows_total <<- n_rows_total + buf_n
    buf <<- NULL
    buf_n <<- 0L
    invisible(TRUE)
  }

  append_rows <- function(df) {
    # df must have the long columns; convert to Arrow Table with stable schema
    if (is.null(df) || nrow(df) < 1) return(invisible(FALSE))
    tbl <- arrow::Table$create(df, schema = long_schema)

    if (is.null(buf)) {
      buf <<- tbl
      buf_n <<- tbl$num_rows
    } else {
      # concat tables
      buf <<- arrow::Table$ConcatTables(list(buf, tbl))
      buf_n <<- buf$num_rows
    }

    if (buf_n >= as.integer(part_rows)) flush_part()
    invisible(TRUE)
  }

  # If no passing models
  if (nrow(passed) < 1) {
    idx <- list(
      schemaVersion = 2,
      generatedAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
      counts = list(
        nModelsTotal = nrow(ms),
        nModelsPassed = 0L,
        nModelsWithTrainingShap = 0L,
        nTrainingSamples = 0L,
        nRows = 0L,
        nParts = 0L
      ),
      artifacts = list(
        training_shap_long_dir = basename(out_dir),
        training_shap_index = basename(out_index_path)
      ),
      models = list()
    )
    jsonlite::write_json(idx, out_index_path, auto_unbox = TRUE, pretty = TRUE)
    if (isTRUE(verbose)) message(glue("[powerup] aggregate_training_shap: no passing models; wrote empty index {out_index_path}"))
    return(invisible(list(ok = TRUE, outputs = list(dir = out_dir, index = out_index_path), counts = idx$counts)))
  }

  # ---- Process each passing model deterministically ----
  for (i in seq_len(nrow(passed))) {
    mk <- passed$modelKey[[i]]
    jobId_i <- passed$jobId[[i]]
    pert <- passed$perturbation[[i]]
    shap_path <- resolve_shap_train_path(passed$shapTrainPath[[i]], mk)

    idx_rec <- list(
      modelKey = mk,
      perturbation = pert,
      passed = TRUE,
      hasTrainingShap = FALSE,
      shapTrainPath = shap_path
    )

    if (!file.exists(shap_path)) {
      if (isTRUE(verbose)) message(glue("[powerup] aggregate_training_shap: missing shap_train.parquet modelKey={mk} path={shap_path}"))
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    # Read parquet as Arrow Table (columnar)
    tab <- tryCatch(
      arrow::read_parquet(shap_path, as_data_frame = FALSE),
      error = function(e) e
    )
    if (inherits(tab, "error") || is.null(tab) || tab$num_rows < 1) {
      if (isTRUE(verbose)) message(glue("[powerup] aggregate_training_shap: empty/invalid shap_train.parquet modelKey={mk} path={shap_path}"))
      if (inherits(tab, "error")) idx_rec$error <- list(type = class(tab)[1], message = conditionMessage(tab))
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    cols <- tab$ColumnNames
    if (length(cols) < 2) {
      idx_rec$error <- list(type = "SchemaError", message = "shap_train.parquet must have cell_line + >=1 feature columns")
      index_models[[length(index_models) + 1L]] <- idx_rec
      next
    }

    # Determine cell_line column
    cell_col <- if ("cell_line" %in% cols) {
      "cell_line"
    } else if ("cellLine" %in% cols) {
      "cellLine"
    } else {
      cols[[1]]
    }

    feature_cols <- setdiff(cols, cell_col)

    # Deterministic feature order: (Intercept) first then alphabetical
    feats <- sort(feature_cols)
    if ("(Intercept)" %in% feats) feats <- c("(Intercept)", setdiff(feats, "(Intercept)"))

    # Exclude intercept rows from output (align with worker behavior)
    feats_out <- setdiff(feats, "(Intercept)")

    idx_rec$trainingShapSchema <- list(
      format = "wide_parquet",
      cell_col = cell_col,
      n_features = length(feats),
      feature_preview = head(feats, 25),
      columns = cols
    )

    # Convert to data.frame in manageable chunks using batches
    # NOTE: read_parquet returns a Table; to_batches() yields record batches.
    batches <- tab$to_batches()

    any_rows <- FALSE

    for (b in batches) {
      df <- as.data.frame(b)

      # Normalize cell_line values
      if (!(cell_col %in% names(df))) next
      cl_vec <- as.character(df[[cell_col]])
      cl_vec <- trimws(cl_vec)
      keep_cl <- !is.na(cl_vec) & nzchar(cl_vec)
      if (!any(keep_cl)) next

      # Track training samples
      for (cl in cl_vec[keep_cl]) train_samples_seen[[cl]] <- TRUE

      # Melt to long using data.table for speed
      dt <- data.table::as.data.table(df)
      data.table::setnames(dt, cell_col, "cell_line")

      # Keep only rows with valid cell_line
      dt <- dt[keep_cl]

      # Melt only feature columns we will output
      # (this avoids intercept and avoids melting unknown extra cols)
      melt_cols <- intersect(feats_out, names(dt))
      if (length(melt_cols) < 1) next

      long_dt <- data.table::melt(
        dt,
        id.vars = "cell_line",
        measure.vars = melt_cols,
        variable.name = "feature",
        value.name = "shap_value",
        variable.factor = FALSE
      )

      # Coerce + filter: keep only abs(shap_value) > 0
      suppressWarnings(long_dt[, shap_value := as.numeric(shap_value)])
      long_dt <- long_dt[!is.na(shap_value) & abs(shap_value) > 0]

      if (nrow(long_dt) < 1) next

      out_df <- data.frame(
        jobId = rep(jobId_i, nrow(long_dt)),
        modelKey = rep(mk, nrow(long_dt)),
        perturbation = rep(pert, nrow(long_dt)),
        cell_line = as.character(long_dt$cell_line),
        feature = as.character(long_dt$feature),
        shap_value = as.numeric(long_dt$shap_value),
        stringsAsFactors = FALSE
      )

      append_rows(out_df)
      any_rows <- TRUE
    }

    if (any_rows) {
      idx_rec$hasTrainingShap <- TRUE
      n_models_with_training_shap <- n_models_with_training_shap + 1L
    }

    index_models[[length(index_models) + 1L]] <- idx_rec

    if (isTRUE(verbose)) message(glue("[powerup] aggregate_training_shap: processed modelKey={mk}"))
  }

  # Flush remaining buffered rows
  flush_part()

  n_train_samples <- length(ls(train_samples_seen, all.names = TRUE))
  n_parts <- length(list.files(out_dir, pattern = "^part-[0-9]+\\.parquet$"))

  idx <- list(
    schemaVersion = 2,
    generatedAt = format(Sys.time(), tz = "UTC", usetz = TRUE),
    counts = list(
      nModelsTotal = nrow(ms),
      nModelsPassed = nrow(passed),
      nModelsWithTrainingShap = n_models_with_training_shap,
      nTrainingSamples = n_train_samples,
      nRows = n_rows_total,
      nParts = n_parts
    ),
    artifacts = list(
      training_shap_long_dir = basename(out_dir),
      training_shap_index = basename(out_index_path)
    ),
    models = index_models
  )

  jsonlite::write_json(idx, out_index_path, auto_unbox = TRUE, pretty = TRUE)

  if (isTRUE(verbose)) {
    message(glue("[powerup] aggregate_training_shap: wrote {out_dir} parts={n_parts} rows={n_rows_total} models_passed={nrow(passed)} models_with_shap={n_models_with_training_shap} training_samples={n_train_samples}"))
    message(glue("[powerup] aggregate_training_shap: wrote {out_index_path} models_indexed={length(index_models)}"))
  }

  invisible(list(
    ok = TRUE,
    outputs = list(dir = out_dir, index = out_index_path),
    counts = idx$counts
  ))
}

