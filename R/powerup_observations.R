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

.pu_obs_safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

.pu_obs_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) < 1) return(NA_character_)
  x[[1]]
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

.pu_obs_normalize_perturbation <- function(x, response_set) {
  strip_ko_prefix <- identical(tolower(response_set), "crispr")
  .pu_canonicalize_perturbation_id(x, strip_ko_prefix = strip_ko_prefix)
}

.pu_obs_load_predictions_required <- function(predictions_path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("[powerup][OBSERVATIONS] arrow R package is required to read parquet predictions.")
  }

  pred_tbl <- arrow::read_parquet(predictions_path) %>% tibble::as_tibble()

  required_cols <- c("cell_line", "perturbation", "modelKey")
  missing_cols <- setdiff(required_cols, colnames(pred_tbl))
  if (length(missing_cols) > 0) {
    stop(glue(
      "[powerup][OBSERVATIONS] predictions parquet missing required columns: ",
      "{paste(missing_cols, collapse=', ')}"
    ))
  }

  # Finalize currently writes pred_mean, but be robust if an older artifact only has pred.
  if (!("pred_mean" %in% colnames(pred_tbl))) {
    if ("pred" %in% colnames(pred_tbl)) {
      pred_tbl <- pred_tbl %>%
        mutate(pred_mean = suppressWarnings(as.numeric(.data$pred)))
    } else {
      stop(
        "[powerup][OBSERVATIONS] predictions parquet must contain pred_mean or pred"
      )
    }
  }

  pred_tbl
}

.pu_obs_extract_target_level <- function(obs_long, response_set) {
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
      normalizedPerturbation = dplyr::if_else(
        is.na(.data$normalizedPerturbation) | !nzchar(trimws(.data$normalizedPerturbation)),
        .pu_obs_normalize_perturbation(.data$perturbation, response_set = response_set),
        as.character(.data$normalizedPerturbation)
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

    # Keep a preferred observation channel for scatter/defaults.
    primary_observation_type <- if (!is.na(target_z)) "target_z" else if (!is.na(mean_target_lfc)) "mean_target_lfc" else .pu_obs_first_nonempty(g$observationType)
    primary_observation_value <- if (!is.na(target_z)) target_z else if (!is.na(mean_target_lfc)) mean_target_lfc else {
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

.pu_obs_resolve_controls <- function(target_tbl, response_set) {
  pkg <- .pu_obs_load_packaged_controls()
  meta <- .pu_obs_load_metadata_controls(target_tbl)

  pos <- unique(c(
    .pu_obs_normalize_perturbation(pkg$positive, response_set),
    meta$positive
  ))
  neg <- unique(c(
    .pu_obs_normalize_perturbation(pkg$negative, response_set),
    meta$negative
  ))

  pos <- sort(unique(pos[!is.na(pos) & nzchar(pos)]))
  neg <- sort(unique(neg[!is.na(neg) & nzchar(neg)]))

  list(
    positive = pos,
    negative = neg,
    source = unique(c(pkg$source, meta$source))
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

.pu_obs_fit_positive_probability <- function(target_tbl) {
  split_tbl <- split(target_tbl, target_tbl$sample, drop = TRUE)

  out <- lapply(split_tbl, function(df) {
    train_df <- df %>%
      filter(!is.na(.data$positive), is.finite(.data$scaled_target_lfc))

    fit_status <- "not_fit"
    fit_message <- NA_character_

    if (nrow(train_df) < 4) {
      fit_status <- "not_fit"
      fit_message <- "Too few control rows to fit logistic model"
      df$positive_probability <- NA_real_
      df$positive_prediction <- NA
      df$positive_probability_model_status <- fit_status
      df$positive_probability_model_message <- fit_message
      return(df)
    }

    if (length(unique(train_df$positive)) < 2) {
      fit_status <- "not_fit"
      fit_message <- "Control rows do not contain both positive and negative classes"
      df$positive_probability <- NA_real_
      df$positive_prediction <- NA
      df$positive_probability_model_status <- fit_status
      df$positive_probability_model_message <- fit_message
      return(df)
    }

    fit <- tryCatch(
      stats::glm(positive ~ scaled_target_lfc, data = train_df, family = "binomial"),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      fit_status <- "failed"
      fit_message <- conditionMessage(fit)
      df$positive_probability <- NA_real_
      df$positive_prediction <- NA
      df$positive_probability_model_status <- fit_status
      df$positive_probability_model_message <- fit_message
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
    df
  })

  bind_rows(out)
}

.pu_obs_build_joined_tbl <- function(target_tbl, pred_tbl, response_set) {
  pred_tbl %>%
    mutate(
      sample = as.character(.data$cell_line),
      normalizedPerturbation = .pu_obs_normalize_perturbation(
        .data$perturbation,
        response_set = response_set
      )
    ) %>%
    left_join(
      target_tbl,
      by = c("sample", "normalizedPerturbation"),
      suffix = c("_pred", "_obs")
    ) %>%
    mutate(
      perturbation_display = dplyr::coalesce(
        .data$perturbation_obs,
        .data$perturbation_pred,
        .data$normalizedPerturbation
      ),
      prediction_value = suppressWarnings(as.numeric(.data$pred_mean))
    )
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

  samples <- target_tbl %>%
    filter(!is.na(.data$sample)) %>%
    distinct(.data$sample) %>%
    arrange(.data$sample) %>%
    transmute(sampleId = .data$sample)

  observation_types <- {
    raw_types <- schema_obj$distinctObservationTypes
    if (is.null(raw_types)) raw_types <- character(0)
    tibble::tibble(observationType = sort(unique(as.character(raw_types))))
  }

  .pu_obs_write_json(
    file.path(out_observations_dir, "samples.json"),
    lapply(seq_len(nrow(samples)), function(i) {
      list(sampleId = as.character(samples$sampleId[[i]]))
    })
  )

  .pu_obs_write_json(
    file.path(out_observations_dir, "observation_types.json"),
    lapply(seq_len(nrow(observation_types)), function(i) {
      list(observationType = as.character(observation_types$observationType[[i]]))
    })
  )

  
  preview_candidates <- joined_tbl %>%
    filter(
      !is.na(.data$primaryObservationType),
      !is.na(.data$primaryObservationValue),
      !is.na(.data$prediction_value)
    ) %>%
    arrange(
      .data$sample,
      .data$primaryObservationType,
      .data$perturbation_display
    )

  scatter_preview <- list(
    schemaVersion = 2,
    generatedAt = as.character(Sys.time()),
    previewSampleId = if (nrow(preview_candidates) > 0) as.character(preview_candidates$sample[[1]]) else NULL,
    previewObservationType = if (nrow(preview_candidates) > 0) as.character(preview_candidates$primaryObservationType[[1]]) else NULL,
    pointCount = min(1000L, nrow(preview_candidates)),
    points = lapply(seq_len(min(1000L, nrow(preview_candidates))), function(i) {
      row <- preview_candidates[i, , drop = FALSE]
      list(
        sampleId = as.character(row$sample[[1]]),
        observationType = as.character(row$primaryObservationType[[1]]),
        perturbation = as.character(row$perturbation_display[[1]]),
        modelKey = as.character(row$modelKey[[1]]),
        x = as.numeric(row$prediction_value[[1]]),
        y = as.numeric(row$primaryObservationValue[[1]]),
        positiveProbability = as.numeric(row$positive_probability[[1]]),
        scaledTargetLfc = as.numeric(row$scaled_target_lfc[[1]])
      )
    })
  )

  .pu_obs_write_json(file.path(out_observations_dir, "manifest.json"), manifest)
  .pu_obs_write_json(file.path(out_observations_dir, "summary.json"), summary)
  .pu_obs_write_json(file.path(out_observations_dir, "scatter_preview.json"), scatter_preview)

  invisible(TRUE)
}

#' Process a submitted observations run.
#'
#' @export
powerup_process_observations <- function(
  cleaned_observations_path,
  schema_path,
  predictions_path,
  out_observations_dir,
  job_id,
  observation_run_id,
  response_set,
  data_version,
  seed = 1L
) {
  powerup_dir_create(out_observations_dir)

  .pu_assert_file_exists(cleaned_observations_path, "cleaned_observations_path")
  .pu_assert_file_exists(schema_path, "schema_path")
  .pu_assert_file_exists(predictions_path, "predictions_path")

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
  pred_tbl <- .pu_obs_load_predictions_required(predictions_path)

  target_tbl <- .pu_obs_extract_target_level(obs_long, response_set = response_set)
  controls <- .pu_obs_resolve_controls(target_tbl, response_set = response_set)

  target_tbl <- target_tbl %>%
    .pu_obs_add_control_annotations(controls = controls) %>%
    .pu_obs_add_scaled_target_lfc()

  positive_tbl <- .pu_obs_fit_positive_probability(target_tbl)

  joined_tbl <- .pu_obs_build_joined_tbl(
    target_tbl = positive_tbl,
    pred_tbl = pred_tbl,
    response_set = response_set
  )

  summary <- list(
    schemaVersion = 2,
    jobId = job_id,
    observationRunId = observation_run_id,
    ok = TRUE,
    mode = "OBSERVATIONS",
    generatedAt = as.character(Sys.time()),
    counts = list(
      nCanonicalRows = nrow(obs_long),
      nTargetLevelRows = nrow(target_tbl),
      nPositiveProbabilityRows = nrow(positive_tbl),
      nJoinedRows = nrow(joined_tbl),
      nSamples = dplyr::n_distinct(target_tbl$sample),
      nObservationTypes = length(unique(obs_long$observationType))
    ),
    controls = list(
      positiveCount = length(controls$positive),
      negativeCount = length(controls$negative),
      source = controls$source
    )
  )

  manifest <- list(
    schemaVersion = 2,
    jobId = job_id,
    observationRunId = observation_run_id,
    mode = "OBSERVATIONS",
    status = "SUCCEEDED",
    responseSet = response_set,
    dataVersion = data_version,
    generatedAt = as.character(Sys.time()),
    inputs = list(
      cleanedObservationsPath = basename(cleaned_observations_path),
      schemaPath = basename(schema_path),
      predictionsPath = basename(predictions_path)
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
    "OBSERVATIONS done target_rows={nrow(target_tbl)} joined_rows={nrow(joined_tbl)}"
  ))

  invisible(TRUE)
}