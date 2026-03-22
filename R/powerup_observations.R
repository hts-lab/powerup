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

  # For cutoff-adjustable posterior work, observations should depend on the
  # cutoff-free Gaussian prior ingredients, not on a precomputed single-cutoff
  # prob_target_event. We therefore require pred_mean and at least one variance
  # representation that can be used later by the frontend/API.
  if (!("pred_mean" %in% colnames(pred_tbl)) && !("pred" %in% colnames(pred_tbl))) {
    stop(
      "[powerup][OBSERVATIONS] predictions parquet must contain pred_mean or pred"
    )
  }

  has_sd_like <- any(c("pred_sd", "pred_var", "pred_log_var") %in% colnames(pred_tbl))
  if (!has_sd_like) {
    stop(
      "[powerup][OBSERVATIONS] predictions parquet must contain at least one of pred_sd, pred_var, or pred_log_var"
    )
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

  # Final type normalization for downstream joins / JSON export
  pred_tbl <- pred_tbl %>%
    mutate(
      pred_mean = suppressWarnings(as.numeric(.data$pred_mean)),
      pred_sd = suppressWarnings(as.numeric(.data$pred_sd)),
      pred_var = suppressWarnings(as.numeric(.data$pred_var)),
      pred_log_var = suppressWarnings(as.numeric(.data$pred_log_var)),
      pred_pi_lower_95 = suppressWarnings(as.numeric(.data$pred_pi_lower_95)),
      pred_pi_upper_95 = suppressWarnings(as.numeric(.data$pred_pi_upper_95)),
      response_cutoff = suppressWarnings(as.numeric(.data$response_cutoff)),
      decreasing = as.logical(.data$decreasing)
    )

  pred_tbl

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
        obsPosteriorEqualPrior = json_num_or_null(row$obs_posterior_equal_prior[[1]])
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
  predictions_path,
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
  .pu_assert_file_exists(predictions_path, "predictions_path")

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
  pred_tbl <- .pu_obs_load_predictions_required(predictions_path)

  target_tbl <- .pu_obs_extract_target_level(
    obs_long,
    response_set = response_set,
    perturbation_tags = perturbation_tags
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
      nObservedSamples = length(observed_samples),
      nPredictedSamples = length(predicted_samples),
      nOverlappingSamples = length(overlapping_samples),
      nObservedPerturbations = length(observed_perts),
      nPredictedPerturbations = length(predicted_perts),
      nOverlappingPerturbations = length(overlapping_perts),
      nObservationTypes = length(unique(obs_long$observationType))
    ),
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
      mode = "not_precomputed_backend",
      reason = "frontend cutoff is user-adjustable, so cutoff-specific prior/posterior must be computed dynamically",
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
    generatedAt = as.character(Sys.time()),
    inputs = list(
      cleanedObservationsPath = basename(cleaned_observations_path),
      schemaPath = basename(schema_path),
      positiveControlsPath = if (!is.null(positive_controls_path) && nzchar(trimws(positive_controls_path))) basename(positive_controls_path) else NULL,
      negativeControlsPath = if (!is.null(negative_controls_path) && nzchar(trimws(negative_controls_path))) basename(negative_controls_path) else NULL,
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
    "OBSERVATIONS done target_rows={nrow(target_tbl)} ",
    "prediction_rows={nrow(pred_tbl)} ",
    "matched_rows={nrow(matched_rows_only)}"
  ))

  invisible(TRUE)
}