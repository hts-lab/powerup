# R/powerup_observations.R

# Internal helpers for PowerUp observation processing live in this file.

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tibble)
  library(glue)
})

.pu_obs_read_json_file <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!nzchar(txt)) return(list())
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

.pu_obs_read_control_file <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(trimws(path))) return(character(0))
  if (!file.exists(path)) return(character(0))

  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!nzchar(trimws(txt))) return(character(0))

  vals <- unlist(strsplit(txt, "[\r\n,;]+", perl = TRUE), use.names = FALSE)
  vals <- trimws(as.character(vals))
  vals <- vals[nzchar(vals)]
  unique(vals)
}


.pu_obs_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) < 1) return(NA_character_)
  x[[1]]
}

.pu_obs_make_clean_name <- function(x) {
  if (length(x) == 0) return(character(0))

  x <- as.character(x)
  x[is.na(x)] <- ""

  x <- trimws(x)
  x <- tolower(x)

  # Replace any run of non-alphanumeric characters with underscore
  x <- gsub("[^a-z0-9]+", "_", x, perl = TRUE)

  # Collapse repeated underscores
  x <- gsub("_+", "_", x, perl = TRUE)

  # Trim leading/trailing underscores
  x <- gsub("^_+", "", x, perl = TRUE)
  x <- gsub("_+$", "", x, perl = TRUE)

  x
}

.pu_obs_parse_meta_row <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(list())
  tryCatch(
    jsonlite::fromJSON(x, simplifyVector = FALSE),
    error = function(e) list(raw = x)
  )
}

.pu_obs_meta_scalar_chr <- function(meta, keys) {
  for (k in keys) {
    v <- meta[[k]]
    if (is.null(v)) next
    if (length(v) < 1) next
    if (is.list(v)) next
    s <- trimws(as.character(v[[1]]))
    if (nzchar(s)) return(s)
  }
  NA_character_
}

.pu_obs_meta_scalar_num <- function(meta, keys) {
  for (k in keys) {
    v <- meta[[k]]
    if (is.null(v)) next
    if (length(v) < 1) next
    if (is.list(v)) next
    n <- suppressWarnings(as.numeric(v[[1]]))
    if (!is.na(n)) return(n)
  }
  NA_real_
}


.pu_obs_strip_known_prefix <- function(x, tags = "ko_") {
  if (length(x) == 0) return(character(0))

  out <- as.character(x)
  tags <- as.character(tags)
  tags <- trimws(tags)
  tags <- unique(tags[nzchar(tags)])

  if (length(tags) < 1) {
    tags <- "ko_"
  }

  for (tag in tags) {
    pattern <- paste0("^", gsub("([][{}()+*^$|\\\\.?-])", "\\\\\\1", tag), collapse = "")
    out <- gsub(pattern, "", out, perl = TRUE, ignore.case = TRUE)
  }

  out
}

.pu_obs_normalize_perturbation <- function(x, response_set, perturbation_tags = "ko_") {
  if (length(x) == 0) return(character(0))

  cleaned <- .pu_obs_make_clean_name(x)

  if (identical(tolower(response_set), "crispr")) {
    tags <- as.character(perturbation_tags)
    tags <- trimws(tags)
    tags <- unique(tags[nzchar(tags)])
    if (length(tags) < 1) {
      tags <- "ko_"
    }

    canonical_tag <- tags[[1]]
    cleaned <- .pu_obs_strip_known_prefix(cleaned, tags = tags)
    return(ifelse(nzchar(cleaned), paste0(canonical_tag, cleaned), cleaned))
  }

  cleaned
}


.pu_obs_normalize_sample_ids <- function(samples) {
  samples <- as.character(samples)
  samples <- trimws(samples)
  samples <- samples[!is.na(samples) & nzchar(samples)]
  unique(samples)
}

.pu_obs_finalize_prediction_tbl <- function(pred_tbl, source_label = "predictions") {
  required_cols <- c("cell_line", "perturbation", "modelKey")
  missing_cols <- setdiff(required_cols, colnames(pred_tbl))
  if (length(missing_cols) > 0) {
    stop(glue(
      "[powerup][OBSERVATIONS] {source_label} missing required columns: ",
      "{paste(missing_cols, collapse=', ')}"
    ))
  }

  if (!("pred_mean" %in% colnames(pred_tbl)) && !("pred" %in% colnames(pred_tbl))) {
    stop(
      glue(
        "[powerup][OBSERVATIONS] {source_label} must contain pred_mean or pred"
      )
    )
  }

  has_sd_like <- any(c("pred_sd", "pred_var", "pred_log_var") %in% colnames(pred_tbl))
  if (!has_sd_like) {
    stop(
      glue(
        "[powerup][OBSERVATIONS] {source_label} must contain at least one of pred_sd, pred_var, or pred_log_var"
      )
    )
  }

  pred_tbl <- pred_tbl %>%
    mutate(
      cell_line = as.character(.data$cell_line),
      perturbation = as.character(.data$perturbation),
      modelKey = as.character(.data$modelKey)
    )

  if (!("pred_mean" %in% colnames(pred_tbl))) {
    if ("pred" %in% colnames(pred_tbl)) {
      pred_tbl <- pred_tbl %>%
        mutate(pred_mean = suppressWarnings(as.numeric(.data$pred)))
    } else {
      stop(
        glue(
          "[powerup][OBSERVATIONS] {source_label} must contain pred_mean or pred"
        )
      )
    }
  }

  if (!("pred_sd" %in% colnames(pred_tbl))) {
    if ("pred_var" %in% colnames(pred_tbl)) {
      pred_tbl <- pred_tbl %>%
        mutate(pred_sd = sqrt(pmax(suppressWarnings(as.numeric(.data$pred_var)), 0)))
    } else if ("pred_log_var" %in% colnames(pred_tbl)) {
      pred_tbl <- pred_tbl %>%
        mutate(pred_sd = sqrt(exp(suppressWarnings(as.numeric(.data$pred_log_var)))))
    } else {
      pred_tbl <- pred_tbl %>%
        mutate(pred_sd = NA_real_)
    }
  }

  if (!("pred_var" %in% colnames(pred_tbl))) {
    pred_tbl <- pred_tbl %>%
      mutate(pred_var = ifelse(is.finite(.data$pred_sd), .data$pred_sd^2, NA_real_))
  }

  if (!("pred_log_var" %in% colnames(pred_tbl))) {
    pred_tbl <- pred_tbl %>%
      mutate(
        pred_log_var = ifelse(
          is.finite(.data$pred_var) & .data$pred_var > 0,
          log(.data$pred_var),
          NA_real_
        )
      )
  }

  if (!("pred_pi_lower_95" %in% colnames(pred_tbl))) {
    pred_tbl <- pred_tbl %>% mutate(pred_pi_lower_95 = NA_real_)
  }

  if (!("pred_pi_upper_95" %in% colnames(pred_tbl))) {
    pred_tbl <- pred_tbl %>% mutate(pred_pi_upper_95 = NA_real_)
  }

  if (!("decreasing" %in% colnames(pred_tbl))) {
    pred_tbl <- pred_tbl %>% mutate(decreasing = NA)
  }

  if (!("response_cutoff" %in% colnames(pred_tbl))) {
    pred_tbl <- pred_tbl %>% mutate(response_cutoff = NA_real_)
  }

  pred_tbl %>%
    mutate(
      pred_mean = suppressWarnings(as.numeric(.data$pred_mean)),
      pred_sd = suppressWarnings(as.numeric(.data$pred_sd)),
      pred_var = suppressWarnings(as.numeric(.data$pred_var)),
      pred_log_var = suppressWarnings(as.numeric(.data$pred_log_var)),
      pred_pi_lower_95 = suppressWarnings(as.numeric(.data$pred_pi_lower_95)),
      pred_pi_upper_95 = suppressWarnings(as.numeric(.data$pred_pi_upper_95)),
      response_cutoff = suppressWarnings(as.numeric(.data$response_cutoff)),
      decreasing = as.logical(.data$decreasing)
    ) %>%
    tibble::as_tibble()
}

.pu_obs_load_predictions_required <- function(
  predictions_path = NULL,
  predictions_duckdb_path = NULL,
  samples = NULL
) {
  samples <- .pu_obs_normalize_sample_ids(samples)

  duckdb_requested <- !is.null(predictions_duckdb_path) &&
    !is.na(predictions_duckdb_path) &&
    nzchar(trimws(predictions_duckdb_path))

  if (duckdb_requested && !file.exists(predictions_duckdb_path)) {
    stop(glue(
      "[powerup][OBSERVATIONS] predictions_duckdb_path was provided but does not exist: {predictions_duckdb_path}"
    ))
  }

  use_duckdb <- duckdb_requested && file.exists(predictions_duckdb_path)
  

  if (use_duckdb) {
    if (!requireNamespace("DBI", quietly = TRUE)) {
      stop("[powerup][OBSERVATIONS] DBI R package is required to read DuckDB predictions.")
    }
    if (!requireNamespace("duckdb", quietly = TRUE)) {
      stop("[powerup][OBSERVATIONS] duckdb R package is required to read DuckDB predictions.")
    }

    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = predictions_duckdb_path, read_only = TRUE)
    on.exit({
      try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    }, add = TRUE)

    if (!DBI::dbExistsTable(con, "user_predictions_long")) {
      stop(
        "[powerup][OBSERVATIONS] DuckDB predictions artifact is missing table 'user_predictions_long'"
      )
    }

    if (length(samples) > 0) {
      DBI::dbWriteTable(
        con,
        "tmp_observation_samples",
        data.frame(cell_line = samples, stringsAsFactors = FALSE),
        temporary = TRUE,
        overwrite = TRUE
      )

      pred_tbl <- DBI::dbGetQuery(
        con,
        "
        SELECT p.*
        FROM user_predictions_long AS p
        INNER JOIN tmp_observation_samples AS s
          ON p.cell_line = s.cell_line
        "
      ) %>%
        tibble::as_tibble()

      message(glue(
        "[powerup][OBSERVATIONS] loaded predictions from DuckDB subset ",
        "samples={length(samples)} rows={nrow(pred_tbl)} db={predictions_duckdb_path}"
      ))
    } else {
      pred_tbl <- DBI::dbGetQuery(
        con,
        "SELECT * FROM user_predictions_long"
      ) %>%
        tibble::as_tibble()

      message(glue(
        "[powerup][OBSERVATIONS] loaded FULL predictions from DuckDB because no sample filter was supplied ",
        "rows={nrow(pred_tbl)} db={predictions_duckdb_path}"
      ))
    }

    return(.pu_obs_finalize_prediction_tbl(pred_tbl, source_label = "predictions DuckDB table"))
  }

  if (is.null(predictions_path) || is.na(predictions_path) || !nzchar(trimws(predictions_path))) {
    stop("[powerup][OBSERVATIONS] either predictions_duckdb_path or predictions_path must be provided")
  }

  .pu_assert_file_exists(predictions_path, "predictions_path")

  pred_tbl <- readr::read_csv(
    predictions_path,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    tibble::as_tibble()

  if (length(samples) > 0 && "cell_line" %in% colnames(pred_tbl)) {
    pred_tbl <- pred_tbl %>%
      filter(.data$cell_line %in% samples)
  }

  message(glue(
    "[powerup][OBSERVATIONS] loaded predictions from CSV ",
    "samples={length(samples)} rows={nrow(pred_tbl)} path={predictions_path}"
  ))

  .pu_obs_finalize_prediction_tbl(pred_tbl, source_label = "predictions CSV")
}


.pu_obs_extract_target_level <- function(obs_long, response_set, perturbation_tags = "ko_") {
  if (!all(c(
    "sample",
    "perturbation",
    "normalizedPerturbation",
    "observationType",
    "observationValue",
    "lineNumber",
    "observationMetaJson"
  ) %in% colnames(obs_long))) {
    stop("[powerup][OBSERVATIONS] cleaned observations CSV missing required columns")
  }

  obs_long <- obs_long %>%
    mutate(
      sample = as.character(.data$sample),
      perturbation = as.character(.data$perturbation),

      # IMPORTANT:
      # Recompute normalizedPerturbation here from perturbation using the
      # R-package-owned canonicalization logic. Do NOT trust the precomputed
      # normalizedPerturbation column from the cleaned CSV, because earlier API
      # parsing may have used different normalization semantics.
      normalizedPerturbation = .pu_obs_normalize_perturbation(
        .data$perturbation,
        response_set = response_set,
        perturbation_tags = perturbation_tags
      ),

      observationType = as.character(.data$observationType),
      observationValue = suppressWarnings(as.numeric(.data$observationValue)),
      lineNumber = suppressWarnings(as.integer(.data$lineNumber)),
      observationMetaJson = as.character(.data$observationMetaJson)
    )

  obs_long$._meta <- lapply(obs_long$observationMetaJson, .pu_obs_parse_meta_row)

  # Build wide-ish target-level table by grouping sample + normalized perturbation.
  grp <- split(obs_long, paste(obs_long$sample, obs_long$normalizedPerturbation, sep = "||"), drop = TRUE)

  out_rows <- vector("list", length(grp))
  idx <- 0L

  for (g in grp) {
    idx <- idx + 1L

    metas <- g$._meta
    meta_first <- metas[[1]]

    sample_id <- .pu_obs_first_nonempty(g$sample)
    normalized_perturbation <- .pu_obs_first_nonempty(g$normalizedPerturbation)
    perturbation_original <- .pu_obs_first_nonempty(g$perturbation)

    # Prefer direct canonical observation rows first; then metadata fallback.
    mean_target_lfc <- {
      vals <- g$observationValue[g$observationType %in% c("avg_lfc", "mean_target_lfc")]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0) vals[[1]] else .pu_obs_meta_scalar_num(meta_first, c("avgLfc", "avg_lfc", "mean_target_lfc"))
    }

    target_z <- {
      vals <- g$observationValue[g$observationType %in% c("z_scored_avg_lfc", "target_z")]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0) vals[[1]] else .pu_obs_meta_scalar_num(meta_first, c("zScoredAvgLfc", "z_scored_avg_lfc", "target_z"))
    }

    target_z_pvalue <- .pu_obs_meta_scalar_num(meta_first, c("zScoredAvgLfcPValue", "z_scored_avg_lfc_p_value", "target_z_pvalue"))
    target_z_fdr <- .pu_obs_meta_scalar_num(meta_first, c("zScoredAvgLfcFdr", "z_scored_avg_lfc_fdr", "target_z_fdr"))
    n_guides <- .pu_obs_meta_scalar_num(meta_first, c("nGuides", "n_guides"))
    control <- .pu_obs_meta_scalar_chr(meta_first, c("control"))
    category_meta <- .pu_obs_meta_scalar_chr(meta_first, c("category", "Category"))
    gene_symbol <- .pu_obs_meta_scalar_chr(meta_first, c("geneSymbol", "gene_symbol", "Gene Symbol"))

    # Keep a preferred observation channel for scatter/defaults using
    # user-facing/raw-style names so filtering stays aligned.
    primary_observation_type <- if (!is.na(target_z)) {
      "z_scored_avg_lfc"
    } else if (!is.na(mean_target_lfc)) {
      "avg_lfc"
    } else {
      .pu_obs_first_nonempty(g$observationType)
    }

    primary_observation_value <- if (!is.na(target_z)) {
      target_z
    } else if (!is.na(mean_target_lfc)) {
      mean_target_lfc
    } else {
      vals <- g$observationValue[is.finite(g$observationValue)]
      if (length(vals) > 0) vals[[1]] else NA_real_
    }

    out_rows[[idx]] <- tibble::tibble(
      sample = sample_id,
      perturbation = perturbation_original,
      normalizedPerturbation = normalized_perturbation,
      gene_symbol = if (!is.na(gene_symbol) && nzchar(gene_symbol)) gene_symbol else perturbation_original,
      mean_target_lfc = mean_target_lfc,
      target_z = target_z,
      target_z_pvalue = target_z_pvalue,
      target_z_fdr = target_z_fdr,
      n_guides = n_guides,
      control = control,
      category_input = category_meta,
      primaryObservationType = primary_observation_type,
      primaryObservationValue = primary_observation_value
    )
  }

  bind_rows(out_rows)
}

.pu_obs_load_packaged_controls <- function() {
  ns <- asNamespace("powerup")

  out <- list(
    positive = character(0),
    negative = character(0),
    source = character(0)
  )

  if (exists("AchillesCommonEssentialControls", envir = ns, inherits = FALSE)) {
    pos <- get("AchillesCommonEssentialControls", envir = ns)
    out$positive <- unique(trimws(as.character(pos)))
    out$source <- c(out$source, "AchillesCommonEssentialControls")
  }

  if (exists("AchillesNonessentialControls", envir = ns, inherits = FALSE)) {
    neg <- get("AchillesNonessentialControls", envir = ns)
    out$negative <- unique(trimws(as.character(neg)))
    out$source <- c(out$source, "AchillesNonessentialControls")
  }

  out
}

.pu_obs_load_metadata_controls <- function(target_tbl) {
  cat_lower <- tolower(trimws(as.character(target_tbl$category_input)))
  ctl_lower <- tolower(trimws(as.character(target_tbl$control)))

  pos <- unique(target_tbl$normalizedPerturbation[
    cat_lower %in% c("positive", "positive_control", "essential", "common_essential") |
      ctl_lower %in% c("positive", "pos", "essential", "common_essential", "true", "1")
  ])
  neg <- unique(target_tbl$normalizedPerturbation[
    cat_lower %in% c("negative", "negative_control", "nonessential", "non_essential") |
      ctl_lower %in% c("negative", "neg", "nonessential", "non_essential", "false", "0")
  ])

  list(
    positive = sort(unique(pos[!is.na(pos) & nzchar(pos)])),
    negative = sort(unique(neg[!is.na(neg) & nzchar(neg)])),
    source = c("observation_metadata")
  )
}

.pu_obs_resolve_controls <- function(
  target_tbl,
  response_set,
  positive_controls_path = NULL,
  negative_controls_path = NULL,
  perturbation_tags = "ko_"
) {
  file_pos_raw <- .pu_obs_read_control_file(positive_controls_path)
  file_neg_raw <- .pu_obs_read_control_file(negative_controls_path)

  pkg <- .pu_obs_load_packaged_controls()
  meta <- .pu_obs_load_metadata_controls(target_tbl)

  file_pos <- .pu_obs_normalize_perturbation(
    file_pos_raw,
    response_set = response_set,
    perturbation_tags = perturbation_tags
  )
  file_neg <- .pu_obs_normalize_perturbation(
    file_neg_raw,
    response_set = response_set,
    perturbation_tags = perturbation_tags
  )

  pkg_pos <- .pu_obs_normalize_perturbation(
    pkg$positive,
    response_set = response_set,
    perturbation_tags = perturbation_tags
  )
  pkg_neg <- .pu_obs_normalize_perturbation(
    pkg$negative,
    response_set = response_set,
    perturbation_tags = perturbation_tags
  )

  pos <- unique(c(file_pos, pkg_pos, meta$positive))
  neg <- unique(c(file_neg, pkg_neg, meta$negative))

  pos <- sort(unique(pos[!is.na(pos) & nzchar(pos)]))
  neg <- sort(unique(neg[!is.na(neg) & nzchar(neg)]))

  list(
    positive = pos,
    negative = neg,
    source = unique(c(
      if (length(file_pos_raw) > 0) "positive_controls_file" else NULL,
      if (length(file_neg_raw) > 0) "negative_controls_file" else NULL,
      pkg$source,
      meta$source
    )),
    uploaded = list(
      positiveRaw = sort(unique(file_pos_raw)),
      negativeRaw = sort(unique(file_neg_raw)),
      positiveNormalized = pos[pos %in% file_pos],
      negativeNormalized = neg[neg %in% file_neg]
    )
  )
}


.pu_obs_add_control_annotations <- function(target_tbl, controls) {
  target_tbl %>%
    mutate(
      control_class = dplyr::case_when(
        .data$normalizedPerturbation %in% controls$positive ~ "positive_control",
        .data$normalizedPerturbation %in% controls$negative ~ "negative_control",
        TRUE ~ "other"
      ),
      positive_control = .data$control_class == "positive_control",
      negative_control = .data$control_class == "negative_control",
      positive = dplyr::case_when(
        .data$positive_control ~ 1L,
        .data$negative_control ~ 0L,
        TRUE ~ NA_integer_
      )
    )
}

.pu_obs_add_scaled_target_lfc <- function(target_tbl) {
  target_tbl %>%
    group_by(.data$sample) %>%
    mutate(
      neg_lfc_median = {
        vals <- .data$mean_target_lfc[.data$negative_control & is.finite(.data$mean_target_lfc)]
        if (length(vals) > 0) stats::median(vals, na.rm = TRUE) else NA_real_
      },
      pos_lfc_median = {
        vals <- .data$mean_target_lfc[.data$positive_control & is.finite(.data$mean_target_lfc)]
        if (length(vals) > 0) stats::median(vals, na.rm = TRUE) else NA_real_
      },
      scaled_target_lfc = dplyr::case_when(
        is.na(.data$mean_target_lfc) ~ NA_real_,
        is.na(.data$neg_lfc_median) ~ NA_real_,
        is.na(.data$pos_lfc_median) ~ NA_real_,
        (.data$neg_lfc_median - .data$pos_lfc_median) == 0 ~ NA_real_,
        TRUE ~ (.data$mean_target_lfc - .data$neg_lfc_median) / (.data$neg_lfc_median - .data$pos_lfc_median)
      )
    ) %>%
    ungroup()
}

.pu_obs_ensure_optional_observation_columns <- function(target_tbl) {
  if (is.null(target_tbl)) {
    return(target_tbl)
  }

  add_na_col <- function(df, col, template) {
    if (!(col %in% colnames(df))) {
      df[[col]] <- template
    }
    df
  }

  n <- nrow(target_tbl)

  # Columns that may legitimately be absent for canonical-long uploads
  # or no-controls runs, but downstream code should always be able to reference.
  target_tbl <- add_na_col(target_tbl, "scaled_target_lfc", rep(NA_real_, n))

  target_tbl <- add_na_col(target_tbl, "positive_probability", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "positive_prediction", rep(NA, n))
  target_tbl <- add_na_col(target_tbl, "positive_probability_model_status", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "positive_probability_model_message", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "positive_probability_model_feature", rep(NA_character_, n))

  target_tbl <- add_na_col(target_tbl, "obs_reps_positive_probability_mean", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_reps_positive_probability_sd", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_reps_positive_probability_n", rep(NA_integer_, n))
  target_tbl <- add_na_col(target_tbl, "obs_reps_positive_probability_status", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "obs_reps_positive_probability_message", rep(NA_character_, n))

  target_tbl <- add_na_col(target_tbl, "obs_stat_used", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "obs_lik_essential", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_lik_nonessential", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_log_lik_essential", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_log_lik_nonessential", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_bayes_factor", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_log_bayes_factor", rep(NA_real_, n))
  target_tbl <- add_na_col(target_tbl, "obs_posterior_equal_prior", rep(NA_real_, n))

  target_tbl <- add_na_col(target_tbl, "obs_likelihood_model", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "obs_likelihood_status", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "obs_likelihood_message", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "obs_n_pos_controls", rep(NA_integer_, n))
  target_tbl <- add_na_col(target_tbl, "obs_n_neg_controls", rep(NA_integer_, n))
  target_tbl <- add_na_col(target_tbl, "obs_kde_bw", rep(NA_character_, n))
  target_tbl <- add_na_col(target_tbl, "obs_kde_adjust", rep(NA_real_, n))

  # Normalize types in case a prior mutate/join created odd types
  target_tbl <- target_tbl %>%
    mutate(
      scaled_target_lfc = suppressWarnings(as.numeric(.data$scaled_target_lfc)),
      positive_probability = suppressWarnings(as.numeric(.data$positive_probability)),
      positive_prediction = as.logical(.data$positive_prediction),
      positive_probability_model_status = as.character(.data$positive_probability_model_status),
      positive_probability_model_message = as.character(.data$positive_probability_model_message),
      positive_probability_model_feature = as.character(.data$positive_probability_model_feature),

      obs_reps_positive_probability_mean = suppressWarnings(as.numeric(.data$obs_reps_positive_probability_mean)),
      obs_reps_positive_probability_sd = suppressWarnings(as.numeric(.data$obs_reps_positive_probability_sd)),
      obs_reps_positive_probability_n = suppressWarnings(as.integer(.data$obs_reps_positive_probability_n)),
      obs_reps_positive_probability_status = as.character(.data$obs_reps_positive_probability_status),
      obs_reps_positive_probability_message = as.character(.data$obs_reps_positive_probability_message),

      obs_stat_used = as.character(.data$obs_stat_used),
      obs_lik_essential = suppressWarnings(as.numeric(.data$obs_lik_essential)),
      obs_lik_nonessential = suppressWarnings(as.numeric(.data$obs_lik_nonessential)),
      obs_log_lik_essential = suppressWarnings(as.numeric(.data$obs_log_lik_essential)),
      obs_log_lik_nonessential = suppressWarnings(as.numeric(.data$obs_log_lik_nonessential)),
      obs_bayes_factor = suppressWarnings(as.numeric(.data$obs_bayes_factor)),
      obs_log_bayes_factor = suppressWarnings(as.numeric(.data$obs_log_bayes_factor)),
      obs_posterior_equal_prior = suppressWarnings(as.numeric(.data$obs_posterior_equal_prior)),
      obs_likelihood_model = as.character(.data$obs_likelihood_model),
      obs_likelihood_status = as.character(.data$obs_likelihood_status),
      obs_likelihood_message = as.character(.data$obs_likelihood_message),
      obs_n_pos_controls = suppressWarnings(as.integer(.data$obs_n_pos_controls)),
      obs_n_neg_controls = suppressWarnings(as.integer(.data$obs_n_neg_controls)),
      obs_kde_bw = as.character(.data$obs_kde_bw),
      obs_kde_adjust = suppressWarnings(as.numeric(.data$obs_kde_adjust))
    )

  target_tbl
}


.pu_obs_fit_positive_probability <- function(target_tbl) {
  split_tbl <- split(target_tbl, target_tbl$sample, drop = TRUE)

  out <- lapply(split_tbl, function(df) {
    # Historical notebook behavior:
    # fit positive ~ mean_target_lfc, not scaled_target_lfc
    train_df <- df %>%
      filter(!is.na(.data$positive), is.finite(.data$mean_target_lfc))

    fit_status <- "not_fit"
    fit_message <- NA_character_

    if (nrow(train_df) < 4) {
      fit_status <- "not_fit"
      fit_message <- "Too few control rows to fit logistic model"
      df$positive_probability <- NA_real_
      df$positive_prediction <- NA
      df$positive_probability_model_status <- fit_status
      df$positive_probability_model_message <- fit_message
      df$positive_probability_model_feature <- "mean_target_lfc"
      return(df)
    }

    if (length(unique(train_df$positive)) < 2) {
      fit_status <- "not_fit"
      fit_message <- "Control rows do not contain both positive and negative classes"
      df$positive_probability <- NA_real_
      df$positive_prediction <- NA
      df$positive_probability_model_status <- fit_status
      df$positive_probability_model_message <- fit_message
      df$positive_probability_model_feature <- "mean_target_lfc"
      return(df)
    }

    fit <- tryCatch(
      stats::glm(positive ~ mean_target_lfc, data = train_df, family = "binomial"),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      fit_status <- "failed"
      fit_message <- conditionMessage(fit)
      df$positive_probability <- NA_real_
      df$positive_prediction <- NA
      df$positive_probability_model_status <- fit_status
      df$positive_probability_model_message <- fit_message
      df$positive_probability_model_feature <- "mean_target_lfc"
      return(df)
    }

    probs <- tryCatch(
      stats::predict(fit, newdata = df, type = "response"),
      error = function(e) rep(NA_real_, nrow(df))
    )

    df$positive_probability <- as.numeric(probs)
    df$positive_prediction <- ifelse(is.na(df$positive_probability), NA, df$positive_probability >= 0.5)
    df$positive_probability_model_status <- "fit_ok"
    df$positive_probability_model_message <- NA_character_
    df$positive_probability_model_feature <- "mean_target_lfc"
    df
  })

  bind_rows(out)
}


.pu_obs_density_eval_safe <- function(x, obs_values, bw = "nrd0", adjust = 1) {
  x <- suppressWarnings(as.numeric(x))
  obs_values <- suppressWarnings(as.numeric(obs_values))
  obs_values <- obs_values[is.finite(obs_values)]

  out <- rep(NA_real_, length(x))

  if (length(obs_values) < 2) {
    return(out)
  }

  # density() can fail if all values are identical or nearly identical.
  # In that case fall back to a narrow Gaussian around the common value.
  if (stats::sd(obs_values, na.rm = TRUE) <= 1e-12) {
    mu <- stats::median(obs_values, na.rm = TRUE)
    sd_fallback <- 0.05
    out[is.finite(x)] <- stats::dnorm(x[is.finite(x)], mean = mu, sd = sd_fallback)
    return(out)
  }

  dens <- tryCatch(
    stats::density(
      obs_values,
      bw = bw,
      adjust = adjust,
      kernel = "gaussian",
      n = 2048,
      from = min(obs_values, na.rm = TRUE) - 4 * stats::sd(obs_values, na.rm = TRUE),
      to   = max(obs_values, na.rm = TRUE) + 4 * stats::sd(obs_values, na.rm = TRUE)
    ),
    error = function(e) NULL
  )

  if (is.null(dens)) {
    return(out)
  }

  finite_idx <- is.finite(x)
  if (!any(finite_idx)) {
    return(out)
  }

  # Linear interpolation, with rule=2 so tails take boundary density
  out[finite_idx] <- approx(
    x = dens$x,
    y = dens$y,
    xout = x[finite_idx],
    rule = 2,
    ties = "ordered"
  )$y

  out
}

.pu_obs_fit_kde_likelihood_by_sample <- function(
  target_tbl,
  obs_col = "target_z",
  bw = "nrd0",
  adjust = 1,
  density_floor = 1e-12
) {
  split_tbl <- split(target_tbl, target_tbl$sample, drop = TRUE)

  out <- lapply(split_tbl, function(df) {
    obs_values <- suppressWarnings(as.numeric(df[[obs_col]]))

    pos_values <- obs_values[df$positive_control & is.finite(obs_values)]
    neg_values <- obs_values[df$negative_control & is.finite(obs_values)]

    df$obs_stat_used <- obs_col
    df$obs_lik_essential <- NA_real_
    df$obs_lik_nonessential <- NA_real_
    df$obs_log_lik_essential <- NA_real_
    df$obs_log_lik_nonessential <- NA_real_
    df$obs_bayes_factor <- NA_real_
    df$obs_log_bayes_factor <- NA_real_
    df$obs_posterior_equal_prior <- NA_real_

    df$obs_likelihood_model <- "kde"
    df$obs_likelihood_status <- "not_fit"
    df$obs_likelihood_message <- NA_character_

    df$obs_n_pos_controls <- length(pos_values)
    df$obs_n_neg_controls <- length(neg_values)
    df$obs_kde_bw <- as.character(bw)
    df$obs_kde_adjust <- as.numeric(adjust)

    if (length(pos_values) < 2 || length(neg_values) < 2) {
      df$obs_likelihood_status <- "not_fit"
      df$obs_likelihood_message <- "Too few finite control values to fit KDE likelihoods"
      return(df)
    }

    lik_pos <- .pu_obs_density_eval_safe(
      x = obs_values,
      obs_values = pos_values,
      bw = bw,
      adjust = adjust
    )
    lik_neg <- .pu_obs_density_eval_safe(
      x = obs_values,
      obs_values = neg_values,
      bw = bw,
      adjust = adjust
    )

    lik_pos <- pmax(as.numeric(lik_pos), density_floor)
    lik_neg <- pmax(as.numeric(lik_neg), density_floor)

    good <- is.finite(obs_values)

    df$obs_lik_essential[good] <- lik_pos[good]
    df$obs_lik_nonessential[good] <- lik_neg[good]
    df$obs_log_lik_essential[good] <- log(lik_pos[good])
    df$obs_log_lik_nonessential[good] <- log(lik_neg[good])
    df$obs_log_bayes_factor[good] <- log(lik_pos[good]) - log(lik_neg[good])
    df$obs_bayes_factor[good] <- exp(df$obs_log_bayes_factor[good])

    # Equal-prior posterior: P(E=1|z) with prior 0.5
    df$obs_posterior_equal_prior[good] <- lik_pos[good] / (lik_pos[good] + lik_neg[good])

    df$obs_likelihood_status <- "fit_ok"
    df$obs_likelihood_message <- NA_character_

    df
  })

  dplyr::bind_rows(out)
}

.pu_obs_add_continuous_posterior_columns <- function(joined_tbl) {
  if (is.null(joined_tbl)) {
    return(joined_tbl)
  }

  add_na_col <- function(df, col, template) {
    if (!(col %in% colnames(df))) {
      df[[col]] <- template
    }
    df
  }

  n <- nrow(joined_tbl)

  joined_tbl <- add_na_col(joined_tbl, "obs_calibration_model", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_calibration_status", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_calibration_reason", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_calibration_intercept", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_calibration_slope", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_calibration_resid_sd", rep(NA_real_, n))

  joined_tbl <- add_na_col(joined_tbl, "obs_implied_response_mean", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_implied_response_sd", rep(NA_real_, n))

  joined_tbl <- add_na_col(joined_tbl, "posterior_mean", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_var", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_sd", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_pi_lower_95", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_pi_upper_95", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_prob_below_cutoff", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_prob_above_cutoff", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_prob_target_event", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_target_event_definition", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_status", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_reason", rep(NA_character_, n))

  joined_tbl <- add_na_col(joined_tbl, "posterior_observation_source", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_observation_value", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "posterior_observation_sd", rep(NA_real_, n))

  joined_tbl <- add_na_col(joined_tbl, "obs_reps_positive_probability_mean", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_reps_positive_probability_sd", rep(NA_real_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_reps_positive_probability_n", rep(NA_integer_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_reps_positive_probability_status", rep(NA_character_, n))
  joined_tbl <- add_na_col(joined_tbl, "obs_reps_positive_probability_message", rep(NA_character_, n))

  joined_tbl <- joined_tbl %>%
    mutate(
      observationType = as.character(.data$observationType),
      prior_pred_mean = suppressWarnings(as.numeric(.data$prior_pred_mean)),
      prior_pred_sd = suppressWarnings(as.numeric(.data$prior_pred_sd)),
      observationValue = suppressWarnings(as.numeric(.data$observationValue)),
      response_cutoff = suppressWarnings(as.numeric(.data$response_cutoff)),
      decreasing = as.logical(.data$decreasing),

      obs_reps_positive_probability_mean = suppressWarnings(as.numeric(.data$obs_reps_positive_probability_mean)),
      obs_reps_positive_probability_sd = suppressWarnings(as.numeric(.data$obs_reps_positive_probability_sd)),
      obs_reps_positive_probability_n = suppressWarnings(as.integer(.data$obs_reps_positive_probability_n)),
      obs_reps_positive_probability_status = as.character(.data$obs_reps_positive_probability_status),
      obs_reps_positive_probability_message = as.character(.data$obs_reps_positive_probability_message),

      obs_calibration_model = as.character(.data$obs_calibration_model),
      obs_calibration_status = as.character(.data$obs_calibration_status),
      obs_calibration_reason = as.character(.data$obs_calibration_reason),
      obs_calibration_intercept = suppressWarnings(as.numeric(.data$obs_calibration_intercept)),
      obs_calibration_slope = suppressWarnings(as.numeric(.data$obs_calibration_slope)),
      obs_calibration_resid_sd = suppressWarnings(as.numeric(.data$obs_calibration_resid_sd)),

      posterior_status = dplyr::if_else(
        is.na(.data$posterior_status) | !nzchar(trimws(.data$posterior_status)),
        "not_available",
        .data$posterior_status
      ),
      posterior_reason = dplyr::if_else(
        is.na(.data$posterior_reason) | !nzchar(trimws(.data$posterior_reason)),
        "No continuous posterior inputs were available.",
        .data$posterior_reason
      )
    )

  prior_mean <- suppressWarnings(as.numeric(joined_tbl$prior_pred_mean))
  prior_sd <- suppressWarnings(as.numeric(joined_tbl$prior_pred_sd))
  cutoff <- suppressWarnings(as.numeric(joined_tbl$response_cutoff))
  decreasing <- as.logical(joined_tbl$decreasing)

  posterior_prob_below_cutoff_fn <- function(mu, sd, cutoff_vec) {
    stats::pnorm(q = cutoff_vec, mean = mu, sd = sd)
  }

  # ---------------------------------------------------------------------------
  # PATH A: same-scale posterior using obs_reps_* values
  # ---------------------------------------------------------------------------
  same_scale_obs_type <- joined_tbl$observationType %in% c(
    "obs_reps_positive_probability_mean"
  )

  same_scale_valid <- (
    same_scale_obs_type &
      is.finite(prior_mean) &
      is.finite(prior_sd) &
      (prior_sd > 0) &
      is.finite(joined_tbl$obs_reps_positive_probability_mean) &
      is.finite(joined_tbl$obs_reps_positive_probability_sd) &
      (joined_tbl$obs_reps_positive_probability_sd > 0)
  )

  if (any(same_scale_valid)) {
    obs_value_same <- joined_tbl$obs_reps_positive_probability_mean[same_scale_valid]
    obs_sd_same <- joined_tbl$obs_reps_positive_probability_sd[same_scale_valid]

    prior_var_same <- prior_sd[same_scale_valid]^2
    obs_var_same <- obs_sd_same^2

    posterior_var_same <- 1 / ((1 / prior_var_same) + (1 / obs_var_same))
    posterior_sd_same <- sqrt(pmax(posterior_var_same, 0))
    posterior_mean_same <- posterior_var_same * (
      (prior_mean[same_scale_valid] / prior_var_same) +
        (obs_value_same / obs_var_same)
    )

    posterior_pi_lower_95_same <- posterior_mean_same - 1.96 * posterior_sd_same
    posterior_pi_upper_95_same <- posterior_mean_same + 1.96 * posterior_sd_same

    posterior_prob_below_cutoff_same <- rep(NA_real_, sum(same_scale_valid))
    posterior_prob_above_cutoff_same <- rep(NA_real_, sum(same_scale_valid))
    posterior_prob_target_event_same <- rep(NA_real_, sum(same_scale_valid))
    posterior_target_event_definition_same <- rep(NA_character_, sum(same_scale_valid))

    cutoff_valid_same <- is.finite(cutoff[same_scale_valid]) &
      is.finite(posterior_mean_same) &
      is.finite(posterior_sd_same) &
      (posterior_sd_same > 0)

    if (any(cutoff_valid_same)) {
      posterior_prob_below_cutoff_same[cutoff_valid_same] <- posterior_prob_below_cutoff_fn(
        mu = posterior_mean_same[cutoff_valid_same],
        sd = posterior_sd_same[cutoff_valid_same],
        cutoff_vec = cutoff[same_scale_valid][cutoff_valid_same]
      )
      posterior_prob_above_cutoff_same[cutoff_valid_same] <- 1 - posterior_prob_below_cutoff_same[cutoff_valid_same]

      dec_vals <- decreasing[same_scale_valid][cutoff_valid_same]
      posterior_prob_target_event_same[cutoff_valid_same] <- ifelse(
        isTRUE(dec_vals) | (!is.na(dec_vals) & dec_vals),
        posterior_prob_below_cutoff_same[cutoff_valid_same],
        posterior_prob_above_cutoff_same[cutoff_valid_same]
      )

      posterior_target_event_definition_same[cutoff_valid_same] <- ifelse(
        isTRUE(dec_vals) | (!is.na(dec_vals) & dec_vals),
        "P(response <= cutoff | prior, obs_reps score)",
        "P(response >= cutoff | prior, obs_reps score)"
      )
    }

    joined_tbl$obs_implied_response_mean[same_scale_valid] <- obs_value_same
    joined_tbl$obs_implied_response_sd[same_scale_valid] <- obs_sd_same

    joined_tbl$posterior_mean[same_scale_valid] <- posterior_mean_same
    joined_tbl$posterior_var[same_scale_valid] <- posterior_var_same
    joined_tbl$posterior_sd[same_scale_valid] <- posterior_sd_same
    joined_tbl$posterior_pi_lower_95[same_scale_valid] <- posterior_pi_lower_95_same
    joined_tbl$posterior_pi_upper_95[same_scale_valid] <- posterior_pi_upper_95_same
    joined_tbl$posterior_prob_below_cutoff[same_scale_valid] <- posterior_prob_below_cutoff_same
    joined_tbl$posterior_prob_above_cutoff[same_scale_valid] <- posterior_prob_above_cutoff_same
    joined_tbl$posterior_prob_target_event[same_scale_valid] <- posterior_prob_target_event_same
    joined_tbl$posterior_target_event_definition[same_scale_valid] <- posterior_target_event_definition_same

    joined_tbl$posterior_observation_source[same_scale_valid] <- "obs_reps_positive_probability"
    joined_tbl$posterior_observation_value[same_scale_valid] <- obs_value_same
    joined_tbl$posterior_observation_sd[same_scale_valid] <- obs_sd_same

    joined_tbl$posterior_status[same_scale_valid] <- "fit_ok_same_scale_obs_reps"
    joined_tbl$posterior_reason[same_scale_valid] <- NA_character_
  }

  same_scale_bad_prior <- (
    same_scale_obs_type &
      !same_scale_valid &
      (
        !is.finite(prior_mean) |
          !is.finite(prior_sd) |
          !(prior_sd > 0)
      )
  )
  joined_tbl$posterior_status[same_scale_bad_prior] <- "not_available_bad_prior"
  joined_tbl$posterior_reason[same_scale_bad_prior] <- "Missing or invalid prior_pred_mean/prior_pred_sd."

  same_scale_bad_obs <- (
    same_scale_obs_type &
      !same_scale_valid &
      !same_scale_bad_prior &
      (
        !is.finite(joined_tbl$obs_reps_positive_probability_mean) |
          !is.finite(joined_tbl$obs_reps_positive_probability_sd) |
          !(joined_tbl$obs_reps_positive_probability_sd > 0)
      )
  )
  joined_tbl$posterior_status[same_scale_bad_obs] <- "not_available_bad_obs_reps"
  joined_tbl$posterior_reason[same_scale_bad_obs] <- paste(
    "obs_reps posterior requires finite obs_reps_positive_probability_mean",
    "and positive obs_reps_positive_probability_sd."
  )

  # ---------------------------------------------------------------------------
  # PATH B: fallback calibration-based posterior for other cases
  # ---------------------------------------------------------------------------
  calibration_candidate <- !same_scale_valid

  has_calibration_cols <- all(c(
    "obs_calibration_intercept",
    "obs_calibration_slope",
    "obs_calibration_resid_sd"
  ) %in% colnames(joined_tbl))

  if (!has_calibration_cols) {
    return(joined_tbl)
  }

  obs_value <- suppressWarnings(as.numeric(joined_tbl$observationValue))
  a <- suppressWarnings(as.numeric(joined_tbl$obs_calibration_intercept))
  b <- suppressWarnings(as.numeric(joined_tbl$obs_calibration_slope))
  obs_sd <- suppressWarnings(as.numeric(joined_tbl$obs_calibration_resid_sd))

  valid <- (
    calibration_candidate &
      is.finite(prior_mean) &
      is.finite(prior_sd) &
      (prior_sd > 0) &
      is.finite(obs_value) &
      is.finite(a) &
      is.finite(b) &
      (abs(b) > 1e-12) &
      is.finite(obs_sd) &
      (obs_sd > 0)
  )

  if (!any(valid)) {
    return(joined_tbl)
  }

  prior_var <- prior_sd[valid]^2
  obs_var <- obs_sd[valid]^2

  obs_implied_mean <- (obs_value[valid] - a[valid]) / b[valid]
  obs_implied_sd <- obs_sd[valid] / abs(b[valid])

  posterior_var <- 1 / ((1 / prior_var) + ((b[valid]^2) / obs_var))
  posterior_sd <- sqrt(pmax(posterior_var, 0))
  posterior_mean <- posterior_var * (
    (prior_mean[valid] / prior_var) +
      (b[valid] * (obs_value[valid] - a[valid]) / obs_var)
  )

  posterior_pi_lower_95 <- posterior_mean - 1.96 * posterior_sd
  posterior_pi_upper_95 <- posterior_mean + 1.96 * posterior_sd

  posterior_prob_below_cutoff <- rep(NA_real_, sum(valid))
  posterior_prob_above_cutoff <- rep(NA_real_, sum(valid))
  posterior_prob_target_event <- rep(NA_real_, sum(valid))
  posterior_target_event_definition <- rep(NA_character_, sum(valid))

  cutoff_valid <- is.finite(cutoff[valid]) & is.finite(posterior_mean) & is.finite(posterior_sd) & (posterior_sd > 0)
  if (any(cutoff_valid)) {
    posterior_prob_below_cutoff[cutoff_valid] <- posterior_prob_below_cutoff_fn(
      mu = posterior_mean[cutoff_valid],
      sd = posterior_sd[cutoff_valid],
      cutoff_vec = cutoff[valid][cutoff_valid]
    )
    posterior_prob_above_cutoff[cutoff_valid] <- 1 - posterior_prob_below_cutoff[cutoff_valid]

    dec_vals <- decreasing[valid][cutoff_valid]
    posterior_prob_target_event[cutoff_valid] <- ifelse(
      isTRUE(dec_vals) | (!is.na(dec_vals) & dec_vals),
      posterior_prob_below_cutoff[cutoff_valid],
      posterior_prob_above_cutoff[cutoff_valid]
    )

    posterior_target_event_definition[cutoff_valid] <- ifelse(
      isTRUE(dec_vals) | (!is.na(dec_vals) & dec_vals),
      "P(response <= cutoff | prior, observation)",
      "P(response >= cutoff | prior, observation)"
    )
  }

  joined_tbl$obs_implied_response_mean[valid] <- obs_implied_mean
  joined_tbl$obs_implied_response_sd[valid] <- obs_implied_sd

  joined_tbl$posterior_mean[valid] <- posterior_mean
  joined_tbl$posterior_var[valid] <- posterior_var
  joined_tbl$posterior_sd[valid] <- posterior_sd
  joined_tbl$posterior_pi_lower_95[valid] <- posterior_pi_lower_95
  joined_tbl$posterior_pi_upper_95[valid] <- posterior_pi_upper_95
  joined_tbl$posterior_prob_below_cutoff[valid] <- posterior_prob_below_cutoff
  joined_tbl$posterior_prob_above_cutoff[valid] <- posterior_prob_above_cutoff
  joined_tbl$posterior_prob_target_event[valid] <- posterior_prob_target_event
  joined_tbl$posterior_target_event_definition[valid] <- posterior_target_event_definition

  joined_tbl$posterior_observation_source[valid] <- "calibrated_observation"
  joined_tbl$posterior_observation_value[valid] <- obs_value[valid]
  joined_tbl$posterior_observation_sd[valid] <- obs_sd[valid]

  joined_tbl$posterior_status[valid] <- "fit_ok_calibrated"
  joined_tbl$posterior_reason[valid] <- NA_character_

  missing_prior <- calibration_candidate & !valid & (
    !is.finite(prior_mean) |
      !is.finite(prior_sd) |
      !(prior_sd > 0)
  )
  joined_tbl$posterior_status[missing_prior] <- "not_available_bad_prior"
  joined_tbl$posterior_reason[missing_prior] <- "Missing or invalid prior_pred_mean/prior_pred_sd."

  missing_obs <- calibration_candidate & !valid & !missing_prior & !is.finite(obs_value)
  joined_tbl$posterior_status[missing_obs] <- "not_available_missing_observation"
  joined_tbl$posterior_reason[missing_obs] <- "Missing or invalid observationValue."

  bad_calibration <- calibration_candidate & !valid & !missing_prior & !missing_obs & (
    !is.finite(a) |
      !is.finite(b) |
      abs(b) <= 1e-12 |
      !is.finite(obs_sd) |
      !(obs_sd > 0)
  )
  joined_tbl$posterior_status[bad_calibration] <- "not_available_bad_calibration"
  joined_tbl$posterior_reason[bad_calibration] <- paste(
    "Continuous posterior requires finite obs_calibration_intercept,",
    "non-zero obs_calibration_slope, and positive obs_calibration_resid_sd."
  )

  joined_tbl
}


.pu_obs_build_observation_long_tbl <- function(target_tbl) {
  target_tbl <- .pu_obs_ensure_optional_observation_columns(target_tbl)
  out <- list()

  # avg_lfc rows
  avg_tbl <- target_tbl %>%
    filter(!is.na(.data$mean_target_lfc)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "avg_lfc",
      observationValue = as.numeric(.data$mean_target_lfc),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- avg_tbl

  # z_scored_avg_lfc rows
  z_tbl <- target_tbl %>%
    filter(!is.na(.data$target_z)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "z_scored_avg_lfc",
      observationValue = as.numeric(.data$target_z),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- z_tbl

  # positive_probability rows
  prob_tbl <- target_tbl %>%
    filter(!is.na(.data$positive_probability)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "positive_probability",
      observationValue = as.numeric(.data$positive_probability),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- prob_tbl


  # obs_reps_positive_probability_mean rows
  obs_reps_prob_tbl <- target_tbl %>%
    filter(!is.na(.data$obs_reps_positive_probability_mean)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "obs_reps_positive_probability_mean",
      observationValue = as.numeric(.data$obs_reps_positive_probability_mean),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_reps_positive_probability_mean = .data$obs_reps_positive_probability_mean,
      obs_reps_positive_probability_sd = .data$obs_reps_positive_probability_sd,
      obs_reps_positive_probability_n = .data$obs_reps_positive_probability_n,
      obs_reps_positive_probability_status = .data$obs_reps_positive_probability_status,
      obs_reps_positive_probability_message = .data$obs_reps_positive_probability_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- obs_reps_prob_tbl


  # obs_lik_essential rows
  lik_ess_tbl <- target_tbl %>%
    filter(!is.na(.data$obs_lik_essential)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "obs_lik_essential",
      observationValue = as.numeric(.data$obs_lik_essential),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- lik_ess_tbl

  # obs_lik_nonessential rows
  lik_non_tbl <- target_tbl %>%
    filter(!is.na(.data$obs_lik_nonessential)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "obs_lik_nonessential",
      observationValue = as.numeric(.data$obs_lik_nonessential),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- lik_non_tbl

  # obs_bayes_factor rows
  bf_tbl <- target_tbl %>%
    filter(!is.na(.data$obs_bayes_factor)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "obs_bayes_factor",
      observationValue = as.numeric(.data$obs_bayes_factor),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- bf_tbl

  # obs_log_bayes_factor rows
  log_bf_tbl <- target_tbl %>%
    filter(!is.na(.data$obs_log_bayes_factor)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "obs_log_bayes_factor",
      observationValue = as.numeric(.data$obs_log_bayes_factor),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- log_bf_tbl

  # obs_posterior_equal_prior rows
  eq_prior_tbl <- target_tbl %>%
    filter(!is.na(.data$obs_posterior_equal_prior)) %>%
    transmute(
      sample = .data$sample,
      perturbation = .data$perturbation,
      normalizedPerturbation = .data$normalizedPerturbation,
      gene_symbol = .data$gene_symbol,
      observationType = "obs_posterior_equal_prior",
      observationValue = as.numeric(.data$obs_posterior_equal_prior),
      control = .data$control,
      category_input = .data$category_input,
      n_guides = .data$n_guides,
      target_z = .data$target_z,
      target_z_pvalue = .data$target_z_pvalue,
      target_z_fdr = .data$target_z_fdr,
      mean_target_lfc = .data$mean_target_lfc,
      scaled_target_lfc = .data$scaled_target_lfc,
      positive_probability = .data$positive_probability,
      positive_prediction = .data$positive_prediction,
      positive_probability_model_status = .data$positive_probability_model_status,
      positive_probability_model_message = .data$positive_probability_model_message,
      obs_stat_used = .data$obs_stat_used,
      obs_lik_essential = .data$obs_lik_essential,
      obs_lik_nonessential = .data$obs_lik_nonessential,
      obs_log_lik_essential = .data$obs_log_lik_essential,
      obs_log_lik_nonessential = .data$obs_log_lik_nonessential,
      obs_bayes_factor = .data$obs_bayes_factor,
      obs_log_bayes_factor = .data$obs_log_bayes_factor,
      obs_posterior_equal_prior = .data$obs_posterior_equal_prior,
      obs_likelihood_model = .data$obs_likelihood_model,
      obs_likelihood_status = .data$obs_likelihood_status,
      obs_likelihood_message = .data$obs_likelihood_message,
      obs_n_pos_controls = .data$obs_n_pos_controls,
      obs_n_neg_controls = .data$obs_n_neg_controls,
      obs_kde_bw = .data$obs_kde_bw,
      obs_kde_adjust = .data$obs_kde_adjust
    )
  out[[length(out) + 1L]] <- eq_prior_tbl

  dplyr::bind_rows(out)
}


.pu_obs_build_joined_tbl <- function(target_tbl, pred_tbl, response_set, perturbation_tags = "ko_") {
  target_tbl <- .pu_obs_ensure_optional_observation_columns(target_tbl)

  obs_long_tbl <- .pu_obs_build_observation_long_tbl(target_tbl) %>%
    mutate(
      sample = .pu_obs_counts_clean_colname(as.character(.data$sample)),
      sample_base = if ("sample_base" %in% colnames(.)) {
        .pu_obs_counts_clean_colname(as.character(.data$sample_base))
      } else {
        .pu_obs_counts_base_sample(as.character(.data$sample))
      },
      normalizedPerturbation = .pu_obs_normalize_perturbation(
        .data$normalizedPerturbation,
        response_set = response_set,
        perturbation_tags = perturbation_tags
      )
    )

    pred_prepped_tbl <- pred_tbl %>%
      mutate(
        sample = .pu_obs_counts_clean_colname(as.character(.data$cell_line)),
        sample_base = .pu_obs_counts_base_sample(as.character(.data$cell_line)),
        normalizedPerturbation = .pu_obs_normalize_perturbation(
          .data$perturbation,
          response_set = response_set,
          perturbation_tags = perturbation_tags
        )
      )

  exact_join_tbl <- pred_prepped_tbl %>%
    inner_join(
      obs_long_tbl,
      by = c("sample", "normalizedPerturbation"),
      suffix = c("_pred", "_obs")
    ) %>%
    mutate(
      sample_match_type = "exact",
      observation_sample = .data$sample,
      observation_sample_base = .data$sample_base_obs,
      prediction_cell_line = as.character(.data$cell_line),
      prediction_sample = .data$sample,
      prediction_sample_base = .data$sample_base_pred
    )

  obs_unmatched_tbl <- obs_long_tbl %>%
    anti_join(
      exact_join_tbl %>%
        distinct(
          .data$sample,
          .data$normalizedPerturbation,
          .data$observationType
        ),
      by = c("sample", "normalizedPerturbation", "observationType")
    )

  fallback_join_tbl <- pred_prepped_tbl %>%
    inner_join(
      obs_unmatched_tbl,
      by = c("sample_base", "normalizedPerturbation"),
      suffix = c("_pred", "_obs")
    ) %>%
    mutate(
      sample_match_type = "base_fallback",
      observation_sample = .data$sample_obs,
      observation_sample_base = .data$sample_base,
      prediction_cell_line = as.character(.data$cell_line),
      prediction_sample = .data$sample_pred,
      prediction_sample_base = .data$sample_base
    )


  joined_tbl <- bind_rows(exact_join_tbl, fallback_join_tbl) %>%
    mutate(
      sampleId = dplyr::coalesce(
        .data$observation_sample,
        .data$sample,
        .data$sample_obs
      ),
      predictionSampleId = dplyr::coalesce(
        .data$prediction_sample,
        .data$sample,
        .data$sample_pred
      )
    )

  required_after_join <- c(
    "perturbation_obs",
    "perturbation_pred",
    "pred_mean",
    "pred_sd",
    "pred_var",
    "pred_log_var",
    "pred_pi_lower_95",
    "pred_pi_upper_95",
    "decreasing",
    "response_cutoff"
  )
  missing_after_join <- setdiff(required_after_join, colnames(joined_tbl))
  if (length(missing_after_join) > 0) {
    stop(glue(
      "[powerup][OBSERVATIONS] joined table missing expected columns after join: ",
      "{paste(missing_after_join, collapse=', ')}. ",
      "Available columns: {paste(colnames(joined_tbl), collapse=', ')}"
    ))
  }

  joined_tbl %>%
    mutate(
      perturbation_display = dplyr::coalesce(
        .data$perturbation_obs,
        .data$perturbation_pred,
        .data$normalizedPerturbation
      ),
      prediction_value = suppressWarnings(as.numeric(.data$pred_mean)),
      prior_pred_mean = suppressWarnings(as.numeric(.data$pred_mean)),
      prior_pred_sd = suppressWarnings(as.numeric(.data$pred_sd)),
      prior_pred_var = suppressWarnings(as.numeric(.data$pred_var)),
      prior_pred_log_var = suppressWarnings(as.numeric(.data$pred_log_var)),
      prior_pred_pi_lower_95 = suppressWarnings(as.numeric(.data$pred_pi_lower_95)),
      prior_pred_pi_upper_95 = suppressWarnings(as.numeric(.data$pred_pi_upper_95)),
      decreasing = .data$decreasing,
      response_cutoff = suppressWarnings(as.numeric(.data$response_cutoff))
    ) %>%
    .pu_obs_add_continuous_posterior_columns()
}




.pu_obs_write_outputs <- function(
  target_tbl,
  positive_tbl,
  joined_tbl,
  out_observations_dir,
  manifest,
  summary,
  schema_obj
) {
  powerup_dir_create(out_observations_dir)

  target_path <- file.path(out_observations_dir, "target_level_observations.parquet")
  positive_path <- file.path(out_observations_dir, "positive_probability.parquet")
  joined_path <- file.path(out_observations_dir, "joined_predictions_observations.parquet")

  .pu_write_parquet_required(target_tbl, target_path)
  .pu_write_parquet_required(positive_tbl, positive_path)
  .pu_write_parquet_required(joined_tbl, joined_path)

  matched_rows_only <- joined_tbl %>%
    filter(
      !is.na(.data$observationType),
      !is.na(.data$observationValue),
      !is.na(.data$prediction_value)
    )

  samples <- target_tbl %>%
    filter(!is.na(.data$sample), nzchar(.data$sample)) %>%
    distinct(.data$sample) %>%
    arrange(.data$sample) %>%
    transmute(sampleId = as.character(.data$sample))

  obs_long_for_schema <- .pu_obs_build_observation_long_tbl(target_tbl)

  observation_types <- obs_long_for_schema %>%
    filter(!is.na(.data$observationType), nzchar(.data$observationType)) %>%
    distinct(.data$observationType) %>%
    arrange(.data$observationType) %>%
    transmute(observationType = as.character(.data$observationType))

  message(glue(
    "[powerup][OBS_WRITE] samples n={nrow(samples)} values={paste(samples$sampleId, collapse=', ')}"
  ))
  message(glue(
    "[powerup][OBS_WRITE] observation_types n={nrow(observation_types)} values={paste(observation_types$observationType, collapse=', ')}"
  ))


  write_json_file <- function(path, obj) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(
      jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      con = path
    )
  }

  json_num_or_null <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) < 1 || is.na(v) || !is.finite(v)) return(NULL)
    v[[1]]
  }

  json_bool_or_null <- function(x) {
    if (length(x) < 1 || is.na(x)) return(NULL)
    as.logical(x[[1]])
  }

  write_json_file(
    file.path(out_observations_dir, "samples.json"),
    lapply(seq_len(nrow(samples)), function(i) {
      list(sampleId = as.character(samples$sampleId[[i]]))
    })
  )

  write_json_file(
    file.path(out_observations_dir, "observation_types.json"),
    lapply(seq_len(nrow(observation_types)), function(i) {
      list(observationType = as.character(observation_types$observationType[[i]]))
    })
  )

  preview_candidates <- joined_tbl %>%
    filter(
      !is.na(.data$observationType),
      !is.na(.data$observationValue),
      !is.na(.data$prediction_value)
    ) %>%
    arrange(
      .data$sampleId,
      .data$observationType,
      .data$perturbation_display
    )

  preview_n <- min(1000L, nrow(preview_candidates))

  scatter_preview <- list(
    schemaVersion = 3,
    generatedAt = as.character(Sys.time()),
    previewSampleId = if (preview_n > 0) as.character(preview_candidates$sampleId[[1]]) else NULL,
    previewObservationType = if (preview_n > 0) as.character(preview_candidates$observationType[[1]]) else NULL,
    pointCount = preview_n,
    points = lapply(seq_len(preview_n), function(i) {
      row <- preview_candidates[i, , drop = FALSE]
      list(
        sampleId = as.character(row$sampleId[[1]]),
        observationType = as.character(row$observationType[[1]]),
        perturbation = as.character(row$perturbation_display[[1]]),
        modelKey = as.character(row$modelKey[[1]]),
        x = json_num_or_null(row$prediction_value[[1]]),
        y = json_num_or_null(row$observationValue[[1]]),

        # Existing experimental score
        positiveProbability = json_num_or_null(row$positive_probability[[1]]),
        scaledTargetLfc = json_num_or_null(row$scaled_target_lfc[[1]]),

        # Cutoff-free prior ingredients from prediction side
        priorPredMean = json_num_or_null(row$prior_pred_mean[[1]]),
        priorPredSd = json_num_or_null(row$prior_pred_sd[[1]]),
        priorPredVar = json_num_or_null(row$prior_pred_var[[1]]),
        priorPredLogVar = json_num_or_null(row$prior_pred_log_var[[1]]),
        priorPredPiLower95 = json_num_or_null(row$prior_pred_pi_lower_95[[1]]),
        priorPredPiUpper95 = json_num_or_null(row$prior_pred_pi_upper_95[[1]]),
        decreasing = json_bool_or_null(row$decreasing[[1]]),
        responseCutoffFromPrediction = json_num_or_null(row$response_cutoff[[1]]),

        # Observation likelihood terms to combine later with any user-selected cutoff
        obsLikEssential = json_num_or_null(row$obs_lik_essential[[1]]),
        obsLikNonessential = json_num_or_null(row$obs_lik_nonessential[[1]]),
        obsBayesFactor = json_num_or_null(row$obs_bayes_factor[[1]]),
        obsLogBayesFactor = json_num_or_null(row$obs_log_bayes_factor[[1]]),
        obsPosteriorEqualPrior = json_num_or_null(row$obs_posterior_equal_prior[[1]]),

        # Continuous posterior layer (additive; may be NA until calibration is available)
        obsCalibrationModel = if ("obs_calibration_model" %in% colnames(row)) as.character(row$obs_calibration_model[[1]]) else NULL,
        obsCalibrationStatus = if ("obs_calibration_status" %in% colnames(row)) as.character(row$obs_calibration_status[[1]]) else NULL,
        obsCalibrationReason = if ("obs_calibration_reason" %in% colnames(row)) as.character(row$obs_calibration_reason[[1]]) else NULL,
        obsCalibrationIntercept = if ("obs_calibration_intercept" %in% colnames(row)) json_num_or_null(row$obs_calibration_intercept[[1]]) else NULL,
        obsCalibrationSlope = if ("obs_calibration_slope" %in% colnames(row)) json_num_or_null(row$obs_calibration_slope[[1]]) else NULL,
        obsCalibrationResidSd = if ("obs_calibration_resid_sd" %in% colnames(row)) json_num_or_null(row$obs_calibration_resid_sd[[1]]) else NULL,
        obsImpliedResponseMean = if ("obs_implied_response_mean" %in% colnames(row)) json_num_or_null(row$obs_implied_response_mean[[1]]) else NULL,
        obsImpliedResponseSd = if ("obs_implied_response_sd" %in% colnames(row)) json_num_or_null(row$obs_implied_response_sd[[1]]) else NULL,

        posteriorMean = if ("posterior_mean" %in% colnames(row)) json_num_or_null(row$posterior_mean[[1]]) else NULL,
        posteriorVar = if ("posterior_var" %in% colnames(row)) json_num_or_null(row$posterior_var[[1]]) else NULL,
        posteriorSd = if ("posterior_sd" %in% colnames(row)) json_num_or_null(row$posterior_sd[[1]]) else NULL,
        posteriorPiLower95 = if ("posterior_pi_lower_95" %in% colnames(row)) json_num_or_null(row$posterior_pi_lower_95[[1]]) else NULL,
        posteriorPiUpper95 = if ("posterior_pi_upper_95" %in% colnames(row)) json_num_or_null(row$posterior_pi_upper_95[[1]]) else NULL,
        posteriorProbBelowCutoff = if ("posterior_prob_below_cutoff" %in% colnames(row)) json_num_or_null(row$posterior_prob_below_cutoff[[1]]) else NULL,
        posteriorProbAboveCutoff = if ("posterior_prob_above_cutoff" %in% colnames(row)) json_num_or_null(row$posterior_prob_above_cutoff[[1]]) else NULL,
        posteriorProbTargetEvent = if ("posterior_prob_target_event" %in% colnames(row)) json_num_or_null(row$posterior_prob_target_event[[1]]) else NULL,
        posteriorTargetEventDefinition = if ("posterior_target_event_definition" %in% colnames(row)) as.character(row$posterior_target_event_definition[[1]]) else NULL,
        posteriorStatus = if ("posterior_status" %in% colnames(row)) as.character(row$posterior_status[[1]]) else NULL,
        posteriorReason = if ("posterior_reason" %in% colnames(row)) as.character(row$posterior_reason[[1]]) else NULL,
        posteriorObservationSource = if ("posterior_observation_source" %in% colnames(row)) as.character(row$posterior_observation_source[[1]]) else NULL,
        posteriorObservationValue = if ("posterior_observation_value" %in% colnames(row)) json_num_or_null(row$posterior_observation_value[[1]]) else NULL,
        posteriorObservationSd = if ("posterior_observation_sd" %in% colnames(row)) json_num_or_null(row$posterior_observation_sd[[1]]) else NULL,

        obsRepsPositiveProbabilityMean = if ("obs_reps_positive_probability_mean" %in% colnames(row)) json_num_or_null(row$obs_reps_positive_probability_mean[[1]]) else NULL,
        obsRepsPositiveProbabilitySd = if ("obs_reps_positive_probability_sd" %in% colnames(row)) json_num_or_null(row$obs_reps_positive_probability_sd[[1]]) else NULL,
        obsRepsPositiveProbabilityN = if ("obs_reps_positive_probability_n" %in% colnames(row)) json_num_or_null(row$obs_reps_positive_probability_n[[1]]) else NULL,
        obsRepsPositiveProbabilityStatus = if ("obs_reps_positive_probability_status" %in% colnames(row)) as.character(row$obs_reps_positive_probability_status[[1]]) else NULL,
        obsRepsPositiveProbabilityMessage = if ("obs_reps_positive_probability_message" %in% colnames(row)) as.character(row$obs_reps_positive_probability_message[[1]]) else NULL


      )
    })
  )

  write_json_file(file.path(out_observations_dir, "manifest.json"), manifest)
  write_json_file(file.path(out_observations_dir, "summary.json"), summary)
  write_json_file(file.path(out_observations_dir, "scatter_preview.json"), scatter_preview)

  invisible(TRUE)
}


#' Process a submitted observations run.
#'
#' @export
powerup_process_observations <- function(
  cleaned_observations_path,
  schema_path,
  positive_controls_path = NULL,
  negative_controls_path = NULL,
  predictions_path = NULL,
  predictions_duckdb_path = NULL,
  out_observations_dir,
  job_id,
  observation_run_id,
  response_set,
  data_version,
  seed = 1L,
  perturbation_tags = "ko_"
) {
  powerup_dir_create(out_observations_dir)

  .pu_assert_file_exists(cleaned_observations_path, "cleaned_observations_path")
  .pu_assert_file_exists(schema_path, "schema_path")

  if (!is.null(predictions_duckdb_path) && !is.na(predictions_duckdb_path) && nzchar(trimws(predictions_duckdb_path))) {
    .pu_assert_file_exists(predictions_duckdb_path, "predictions_duckdb_path")
  } else if (!is.null(predictions_path) && !is.na(predictions_path) && nzchar(trimws(predictions_path))) {
    .pu_assert_file_exists(predictions_path, "predictions_path")
  } else {
    stop("[powerup][OBSERVATIONS] either predictions_duckdb_path or predictions_path must be provided")
  }

  if (!is.null(positive_controls_path) && !is.na(positive_controls_path) && nzchar(trimws(positive_controls_path))) {
    .pu_assert_file_exists(positive_controls_path, "positive_controls_path")
  }
  if (!is.null(negative_controls_path) && !is.na(negative_controls_path) && nzchar(trimws(negative_controls_path))) {
    .pu_assert_file_exists(negative_controls_path, "negative_controls_path")
  }

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "OBSERVATIONS start response_set={response_set} data_version={data_version} seed={seed}"
  ))

  obs_long <- readr::read_csv(
    cleaned_observations_path,
    show_col_types = FALSE,
    progress = FALSE
  ) %>% tibble::as_tibble()

  schema_obj <- .pu_obs_read_json_file(schema_path)

  target_tbl <- .pu_obs_extract_target_level(
    obs_long,
    response_set = response_set,
    perturbation_tags = perturbation_tags
  )

  observation_samples <- unique(trimws(as.character(target_tbl$sample)))
  observation_samples <- observation_samples[!is.na(observation_samples) & nzchar(observation_samples)]

  pred_tbl <- .pu_obs_load_predictions_required(
    predictions_path = predictions_path,
    predictions_duckdb_path = predictions_duckdb_path,
    samples = observation_samples
  )

  controls <- .pu_obs_resolve_controls(
    target_tbl = target_tbl,
    response_set = response_set,
    positive_controls_path = positive_controls_path,
    negative_controls_path = negative_controls_path,
    perturbation_tags = perturbation_tags
  )

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "controls positive={length(controls$positive)} negative={length(controls$negative)} ",
    "sources={paste(controls$source, collapse='|')}"
  ))

  if (!is.null(controls$uploaded)) {
    message(glue(
      "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
      "uploaded_control_counts raw_positive={length(controls$uploaded$positiveRaw)} ",
      "raw_negative={length(controls$uploaded$negativeRaw)} ",
      "normalized_positive={length(controls$uploaded$positiveNormalized)} ",
      "normalized_negative={length(controls$uploaded$negativeNormalized)}"
    ))
  }

  target_tbl <- target_tbl %>%
    .pu_obs_add_control_annotations(controls = controls) %>%
    .pu_obs_add_scaled_target_lfc() %>%
    .pu_obs_ensure_optional_observation_columns()

  # Keep existing discriminative score
  positive_tbl <- .pu_obs_fit_positive_probability(target_tbl) %>%
    .pu_obs_ensure_optional_observation_columns()


  fit_feature_examples <- positive_tbl %>%
    dplyr::distinct(.data$sample, .data$positive_probability_model_feature, .data$positive_probability_model_status) %>%
    dplyr::arrange(.data$sample)

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "positive_probability_features={paste(unique(fit_feature_examples$positive_probability_model_feature), collapse='|')} ",
    "fit_statuses={paste(unique(fit_feature_examples$positive_probability_model_status), collapse='|')}"
  ))

  # Add generative KDE likelihood model using target_z
  positive_tbl <- .pu_obs_fit_kde_likelihood_by_sample(
    target_tbl = positive_tbl,
    obs_col = "target_z",
    bw = "nrd0",
    adjust = 1,
    density_floor = 1e-12
  ) %>%
    .pu_obs_ensure_optional_observation_columns()

  joined_tbl <- .pu_obs_build_joined_tbl(
    target_tbl = positive_tbl,
    pred_tbl = pred_tbl,
    response_set = response_set,
    perturbation_tags = perturbation_tags
  )

  pred_norm_tbl <- pred_tbl %>%
    mutate(
      normalizedPerturbationObsJoin = .pu_obs_normalize_perturbation(
        .data$perturbation,
        response_set = response_set,
        perturbation_tags = perturbation_tags
      )
    )

  observed_perts <- sort(unique(as.character(
    positive_tbl$normalizedPerturbation[!is.na(positive_tbl$normalizedPerturbation)]
  )))
  predicted_perts <- sort(unique(as.character(
    pred_norm_tbl$normalizedPerturbationObsJoin[!is.na(pred_norm_tbl$normalizedPerturbationObsJoin)]
  )))
  overlapping_perts <- intersect(observed_perts, predicted_perts)

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "observed_perturbations={length(observed_perts)} ",
    "predicted_perturbations={length(predicted_perts)} ",
    "overlapping_perturbations={length(overlapping_perts)}"
  ))

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "observed_perturbation_examples={paste(head(observed_perts, 20), collapse=', ')}"
  ))

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "predicted_perturbation_examples={paste(head(predicted_perts, 20), collapse=', ')}"
  ))

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "overlapping_perturbation_examples={paste(head(overlapping_perts, 20), collapse=', ')}"
  ))

  matched_rows_only <- joined_tbl %>%
    filter(
      !is.na(.data$observationType),
      !is.na(.data$observationValue),
      !is.na(.data$prediction_value)
    )

  nPriorGaussianRows <- sum(
    is.finite(suppressWarnings(as.numeric(joined_tbl$prior_pred_mean))) &
      is.finite(suppressWarnings(as.numeric(joined_tbl$prior_pred_sd))),
    na.rm = TRUE
  )

  nObsLikelihoodRows <- sum(
    is.finite(suppressWarnings(as.numeric(joined_tbl$obs_lik_essential))) &
      is.finite(suppressWarnings(as.numeric(joined_tbl$obs_lik_nonessential))),
    na.rm = TRUE
  )


  nContinuousPosteriorRows <- sum(
    is.finite(suppressWarnings(as.numeric(joined_tbl$posterior_mean))) &
      is.finite(suppressWarnings(as.numeric(joined_tbl$posterior_sd))),
    na.rm = TRUE
  )

observed_samples <- sort(unique(as.character(
  target_tbl$sample[!is.na(target_tbl$sample)]
)))
predicted_samples <- sort(unique(as.character(
  pred_tbl$cell_line[!is.na(pred_tbl$cell_line)]
)))

observed_samples_clean <- sort(unique(.pu_obs_counts_clean_colname(observed_samples)))
predicted_samples_clean <- sort(unique(.pu_obs_counts_clean_colname(predicted_samples)))

observed_sample_bases <- sort(unique(.pu_obs_counts_base_sample(observed_samples)))
predicted_sample_bases <- sort(unique(.pu_obs_counts_base_sample(predicted_samples)))

overlapping_samples <- intersect(observed_samples_clean, predicted_samples_clean)
overlapping_sample_bases <- intersect(observed_sample_bases, predicted_sample_bases)


  summary <- list(
    schemaVersion = 3,
    jobId = job_id,
    observationRunId = observation_run_id,
    ok = TRUE,
    mode = "OBSERVATIONS",
    generatedAt = as.character(Sys.time()),
    counts = list(
      nCanonicalRows = nrow(obs_long),
      nTargetLevelRows = nrow(target_tbl),
      nPositiveProbabilityRows = nrow(positive_tbl),
      nPredictionRows = nrow(pred_tbl),
      nMatchedRows = nrow(matched_rows_only),
      nPriorGaussianRows = nPriorGaussianRows,
      nObsLikelihoodRows = nObsLikelihoodRows,
      nContinuousPosteriorRows = nContinuousPosteriorRows,
      nObservedSamples = length(observed_samples),
      nPredictedSamples = length(predicted_samples),
      nOverlappingSamples = length(overlapping_samples),
      nObservedPerturbations = length(observed_perts),
      nPredictedPerturbations = length(predicted_perts),
      nOverlappingPerturbations = length(overlapping_perts),
      nObservationTypes = length(unique(obs_long$observationType))
    ),
    predictionsSource = if (!is.null(predictions_duckdb_path) && nzchar(trimws(predictions_duckdb_path))) "duckdb" else "csv",
    sampleOverlap = list(
      overlappingExact = head(overlapping_samples, 50),
      overlappingBase = head(overlapping_sample_bases, 50),
      observedOnlyExact = head(setdiff(observed_samples_clean, predicted_samples_clean), 50),
      predictedOnlyExact = head(setdiff(predicted_samples_clean, observed_samples_clean), 50),
      observedBasesOnly = head(setdiff(observed_sample_bases, predicted_sample_bases), 50),
      predictedBasesOnly = head(setdiff(predicted_sample_bases, observed_sample_bases), 50)
    ),
    perturbationOverlap = list(
      overlapping = head(overlapping_perts, 50),
      observedOnly = head(setdiff(observed_perts, predicted_perts), 50),
      predictedOnly = head(setdiff(predicted_perts, observed_perts), 50)
    ),
    controls = list(
      positiveCount = length(controls$positive),
      negativeCount = length(controls$negative),
      source = controls$source,
      uploadedPositiveRawCount = if (!is.null(controls$uploaded)) length(controls$uploaded$positiveRaw) else 0L,
      uploadedNegativeRawCount = if (!is.null(controls$uploaded)) length(controls$uploaded$negativeRaw) else 0L,
      uploadedPositiveNormalizedCount = if (!is.null(controls$uploaded)) length(controls$uploaded$positiveNormalized) else 0L,
      uploadedNegativeNormalizedCount = if (!is.null(controls$uploaded)) length(controls$uploaded$negativeNormalized) else 0L
    ),
    posterior = list(
      mode = "continuous_posterior_columns_added",
      reason = paste(
        "Backend now writes additive continuous-posterior columns using either",
        "a same-scale observation path (for example obs_reps-positive-probability",
        "mean/sd when available) or a continuous calibration path when available.",
        "Current KDE class-likelihood fields remain unchanged and are still written."
      ),
      priorDistribution = list(
        family = "gaussian",
        meanColumn = "prior_pred_mean",
        sdColumn = "prior_pred_sd",
        varColumn = "prior_pred_var",
        logVarColumn = "prior_pred_log_var",
        piLower95Column = "prior_pred_pi_lower_95",
        piUpper95Column = "prior_pred_pi_upper_95",
        decreasingColumn = "decreasing"
      ),
      observationLikelihood = list(
        essentialColumn = "obs_lik_essential",
        nonessentialColumn = "obs_lik_nonessential",
        bayesFactorColumn = "obs_bayes_factor",
        logBayesFactorColumn = "obs_log_bayes_factor",
        equalPriorPosteriorColumn = "obs_posterior_equal_prior",
        observationStatisticUsed = "target_z",
        model = "kde"
      ),
      continuousPosterior = list(
        family = "gaussian",
        statusColumn = "posterior_status",
        reasonColumn = "posterior_reason",
        meanColumn = "posterior_mean",
        varColumn = "posterior_var",
        sdColumn = "posterior_sd",
        piLower95Column = "posterior_pi_lower_95",
        piUpper95Column = "posterior_pi_upper_95",
        probBelowCutoffColumn = "posterior_prob_below_cutoff",
        probAboveCutoffColumn = "posterior_prob_above_cutoff",
        probTargetEventColumn = "posterior_prob_target_event",
        targetEventDefinitionColumn = "posterior_target_event_definition",
        observationCalibrationModelColumn = "obs_calibration_model",
        observationCalibrationStatusColumn = "obs_calibration_status",
        observationCalibrationReasonColumn = "obs_calibration_reason",
        observationCalibrationInterceptColumn = "obs_calibration_intercept",
        observationCalibrationSlopeColumn = "obs_calibration_slope",
        observationCalibrationResidSdColumn = "obs_calibration_resid_sd",
        observationImpliedResponseMeanColumn = "obs_implied_response_mean",
        observationImpliedResponseSdColumn = "obs_implied_response_sd",
        observationSourceColumn = "posterior_observation_source",
        observationValueColumn = "posterior_observation_value",
        observationSdColumn = "posterior_observation_sd",
        obsRepsMeanColumn = "obs_reps_positive_probability_mean",
        obsRepsSdColumn = "obs_reps_positive_probability_sd",
        obsRepsNColumn = "obs_reps_positive_probability_n",
        obsRepsStatusColumn = "obs_reps_positive_probability_status",
        obsRepsMessageColumn = "obs_reps_positive_probability_message"
      )
    )
  )

  manifest <- list(
    schemaVersion = 3,
    jobId = job_id,
    observationRunId = observation_run_id,
    mode = "OBSERVATIONS",
    status = "SUCCEEDED",
    responseSet = response_set,
    dataVersion = data_version,
    predictionsSource = if (!is.null(predictions_duckdb_path) && nzchar(trimws(predictions_duckdb_path))) "duckdb" else "csv",
    generatedAt = as.character(Sys.time()),
    inputs = list(
      cleanedObservationsPath = basename(cleaned_observations_path),
      schemaPath = basename(schema_path),
      positiveControlsPath = if (!is.null(positive_controls_path) && nzchar(trimws(positive_controls_path))) basename(positive_controls_path) else NULL,
      negativeControlsPath = if (!is.null(negative_controls_path) && nzchar(trimws(negative_controls_path))) basename(negative_controls_path) else NULL,
      predictionsPath = if (!is.null(predictions_path) && nzchar(trimws(predictions_path))) basename(predictions_path) else NULL,
      predictionsDuckdbPath = if (!is.null(predictions_duckdb_path) && nzchar(trimws(predictions_duckdb_path))) basename(predictions_duckdb_path) else NULL
    ),
    artifacts = list(
      manifest = "manifest.json",
      summary = "summary.json",
      scatterPreview = "scatter_preview.json",
      targetLevelObservations = "target_level_observations.parquet",
      positiveProbability = "positive_probability.parquet",
      joinedPredictionsObservations = "joined_predictions_observations.parquet",
      samples = "samples.json",
      observationTypes = "observation_types.json"
    )
  )

  .pu_obs_write_outputs(
    target_tbl = target_tbl,
    positive_tbl = positive_tbl,
    joined_tbl = joined_tbl,
    out_observations_dir = out_observations_dir,
    manifest = manifest,
    summary = summary,
    schema_obj = schema_obj
  )

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "OBSERVATIONS done target_rows={nrow(target_tbl)} ",
    "prediction_rows={nrow(pred_tbl)} ",
    "matched_rows={nrow(matched_rows_only)}"
  ))

  invisible(TRUE)
}