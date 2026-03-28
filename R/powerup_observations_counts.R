# R/powerup_observations_counts.R

# Raw-counts CRISPR observations pipeline for PowerUp.
# This path is intended to reproduce the poola-based workflow as closely
# as possible, starting from guide-level raw counts rather than target-level
# canonical long observations.
# Credit to the original poola implementation at https://github.com/broadinstitute/poola

suppressPackageStartupMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(glue)
})

# -----------------------------
# Counts-mode helpers
# -----------------------------

.pu_obs_counts_find_col <- function(colnames_vec, candidates) {
  if (length(colnames_vec) < 1) {
    return(NA_character_)
  }

  raw_names <- as.character(colnames_vec)
  norm_names <- gsub("[^a-z0-9]+", "", tolower(trimws(raw_names)))

  candidate_norm <- gsub("[^a-z0-9]+", "", tolower(trimws(as.character(candidates))))
  candidate_norm <- unique(candidate_norm[nzchar(candidate_norm)])

  hit_idx <- match(candidate_norm, norm_names, nomatch = 0L)
  hit_idx <- hit_idx[hit_idx > 0]

  if (length(hit_idx) < 1) {
    return(NA_character_)
  }

  raw_names[[hit_idx[[1]]]]
}

.pu_obs_counts_guess_delim <- function(path) {
  .pu_assert_file_exists(path, "path")

  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) < 1) {
    stop("[powerup][OBSERVATIONS_COUNTS] input file is empty")
  }

  line <- first_line[[1]]

  n_tab <- lengths(regmatches(line, gregexpr("\t", line, fixed = TRUE)))
  n_comma <- lengths(regmatches(line, gregexpr(",", line, fixed = TRUE)))
  n_semicolon <- lengths(regmatches(line, gregexpr(";", line, fixed = TRUE)))

  n_tab <- ifelse(length(n_tab) < 1 || is.na(n_tab), 0L, as.integer(n_tab[[1]]))
  n_comma <- ifelse(length(n_comma) < 1 || is.na(n_comma), 0L, as.integer(n_comma[[1]]))
  n_semicolon <- ifelse(length(n_semicolon) < 1 || is.na(n_semicolon), 0L, as.integer(n_semicolon[[1]]))

  if (n_tab >= n_comma && n_tab >= n_semicolon && n_tab > 0) {
    return("\t")
  }
  if (n_comma >= n_semicolon && n_comma > 0) {
    return(",")
  }
  if (n_semicolon > 0) {
    return(";")
  }

  stop(
    "[powerup][OBSERVATIONS_COUNTS] could not determine delimiter from first line; ",
    "expected tab, comma, or semicolon delimited text"
  )
}

.pu_obs_counts_clean_colname <- function(x) {
  x <- as.character(x)
  x <- .pu_obs_make_clean_name(x)
  x <- as.character(x)
  x
}

.pu_obs_counts_base_sample <- function(x) {
  x <- .pu_obs_counts_clean_colname(x)
  x <- gsub("_[^_]+$", "", x, perl = TRUE)
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

.pu_obs_counts_step <- function(step_name, expr) {
  message(glue("[powerup][OBSERVATIONS_COUNTS] STEP_START {step_name}"))
  out <- tryCatch(
    force(expr),
    error = function(e) {
      message(glue(
        "[powerup][OBSERVATIONS_COUNTS] STEP_FAIL {step_name} err={conditionMessage(e)}"
      ))
      stop(e)
    }
  )
  message(glue("[powerup][OBSERVATIONS_COUNTS] STEP_DONE {step_name}"))
  out
}


.pu_obs_counts_read_reference_required <- function(reference_path) {
  .pu_assert_file_exists(reference_path, "reference_path")

  delim <- .pu_obs_counts_guess_delim(reference_path)

  ref_tbl <- readr::read_delim(
    file = reference_path,
    delim = delim,
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 10000,
    trim_ws = TRUE
  ) %>%
    tibble::as_tibble()

  if (ncol(ref_tbl) < 2) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] guide reference parsed into only {ncol(ref_tbl)} column(s) ",
      "using delimiter '{delim}'. Header may be malformed."
    ))
  }

  original_colnames <- colnames(ref_tbl)
  colnames(ref_tbl) <- .pu_obs_make_clean_name(colnames(ref_tbl))
  cleaned_colnames <- colnames(ref_tbl)

  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] guide_reference delimiter='{delim}'"
  ))
  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] guide_reference original columns: {paste(original_colnames, collapse=', ')}"
  ))
  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] guide_reference cleaned columns: {paste(cleaned_colnames, collapse=', ')}"
  ))

  barcode_col <- .pu_obs_counts_find_col(
    colnames(ref_tbl),
    c(
      "barcode_sequence",
      "construct_barcode",
      "barcodesequence",
      "constructbarcode",
      "barcode sequence",
      "construct barcode"
    )
  )

  gene_col <- .pu_obs_counts_find_col(
    colnames(ref_tbl),
    c(
      "gene_symbol",
      "gene",
      "genesymbol",
      "gene symbol"
    )
  )

  gene_id_col <- .pu_obs_counts_find_col(
    colnames(ref_tbl),
    c(
      "gene_id",
      "geneid",
      "gene id"
    )
  )

  if (is.na(barcode_col) || !nzchar(barcode_col)) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] reference file missing barcode_sequence/construct_barcode column. ",
      "Available cleaned columns: {paste(colnames(ref_tbl), collapse=', ')}"
    ))
  }

  if (is.na(gene_col) || !nzchar(gene_col)) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] reference file missing gene_symbol/gene column. ",
      "Available cleaned columns: {paste(colnames(ref_tbl), collapse=', ')}"
    ))
  }

  ref_tbl <- ref_tbl %>%
    transmute(
      construct_barcode = as.character(.data[[barcode_col]]),
      gene_symbol = as.character(.data[[gene_col]]),
      gene_id = if (!is.na(gene_id_col) && nzchar(gene_id_col)) {
        as.character(.data[[gene_id_col]])
      } else {
        NA_character_
      }
    ) %>%
    filter(!is.na(.data$construct_barcode), nzchar(.data$construct_barcode)) %>%
    distinct()

  ref_tbl
}


.pu_obs_counts_read_counts_required <- function(counts_path) {
  .pu_assert_file_exists(counts_path, "counts_path")

  delim <- .pu_obs_counts_guess_delim(counts_path)

  counts_tbl <- readr::read_delim(
    file = counts_path,
    delim = delim,
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 10000,
    trim_ws = TRUE
  ) %>%
    tibble::as_tibble()

  if (ncol(counts_tbl) < 2) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] raw counts file parsed into only {ncol(counts_tbl)} column(s) ",
      "using delimiter '{delim}'. Header may be malformed."
    ))
  }

  original_colnames <- colnames(counts_tbl)
  colnames(counts_tbl) <- .pu_obs_make_clean_name(colnames(counts_tbl))
  cleaned_colnames <- colnames(counts_tbl)

  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] raw_counts delimiter='{delim}'"
  ))
  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] raw_counts original columns: {paste(original_colnames, collapse=', ')}"
  ))
  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] raw_counts cleaned columns: {paste(cleaned_colnames, collapse=', ')}"
  ))

  barcode_col <- .pu_obs_counts_find_col(
    colnames(counts_tbl),
    c(
      "construct_barcode",
      "barcode_sequence",
      "constructbarcode",
      "barcodesequence",
      "construct barcode",
      "barcode sequence"
    )
  )

  construct_id_col <- .pu_obs_counts_find_col(
    colnames(counts_tbl),
    c(
      "construct_id",
      "construct_ids",
      "construct_i_ds",
      "constructid",
      "constructids",
      "construct id",
      "construct ids"
    )
  )

  if (is.na(barcode_col) || !nzchar(barcode_col)) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] counts file missing construct_barcode/barcode_sequence column. ",
      "Available cleaned columns: {paste(colnames(counts_tbl), collapse=', ')}"
    ))
  }

  counts_tbl <- counts_tbl %>%
    rename(construct_barcode = all_of(barcode_col))

  if (!is.na(construct_id_col) && nzchar(construct_id_col) && construct_id_col != "construct_id") {
    counts_tbl <- counts_tbl %>%
      rename(construct_id = all_of(construct_id_col))
  } else if (!("construct_id" %in% colnames(counts_tbl))) {
    counts_tbl <- counts_tbl %>%
      mutate(construct_id = NA_character_)
  }

  counts_tbl <- counts_tbl %>%
    mutate(
      construct_barcode = as.character(.data$construct_barcode),
      construct_id = as.character(.data$construct_id)
    )

  counts_tbl
}


.pu_obs_counts_parse_config <- function(config_path = NULL) {
  cfg <- list()

  if (!is.null(config_path) && length(config_path) >= 1) {
    config_path <- as.character(config_path[[1]])
    if (!is.na(config_path)) {
      config_path <- trimws(config_path)
      if (nzchar(config_path)) {
        .pu_assert_file_exists(config_path, "config_path")
        cfg <- .pu_obs_read_json_file(config_path)
      }
    }
  }

  # --- guideNegativeControlPatterns ---
  guide_negative_control_patterns <- cfg[["guideNegativeControlPatterns"]]
  if (is.null(guide_negative_control_patterns) || length(guide_negative_control_patterns) < 1) {
    guide_negative_control_patterns <- c("ONE_INTERGENIC_SITE", "NO_SITE")
  }
  guide_negative_control_patterns <- as.character(
    unlist(guide_negative_control_patterns, use.names = FALSE)
  )
  guide_negative_control_patterns <- guide_negative_control_patterns[
    !is.na(guide_negative_control_patterns) & nzchar(trimws(guide_negative_control_patterns))
  ]
  if (length(guide_negative_control_patterns) < 1) {
    guide_negative_control_patterns <- c("ONE_INTERGENIC_SITE", "NO_SITE")
  }

  # --- guideNegativeControlFuzzy ---
  guide_negative_control_fuzzy <- cfg[["guideNegativeControlFuzzy"]]
  if (is.null(guide_negative_control_fuzzy) || length(guide_negative_control_fuzzy) < 1) {
    guide_negative_control_fuzzy <- rep(TRUE, length(guide_negative_control_patterns))
  }
  guide_negative_control_fuzzy <- as.logical(
    unlist(guide_negative_control_fuzzy, use.names = FALSE)
  )
  if (length(guide_negative_control_fuzzy) < 1) {
    guide_negative_control_fuzzy <- rep(TRUE, length(guide_negative_control_patterns))
  }
  if (length(guide_negative_control_fuzzy) < length(guide_negative_control_patterns)) {
    guide_negative_control_fuzzy <- c(
      guide_negative_control_fuzzy,
      rep(TRUE, length(guide_negative_control_patterns) - length(guide_negative_control_fuzzy))
    )
  }
  guide_negative_control_fuzzy <- guide_negative_control_fuzzy[seq_len(length(guide_negative_control_patterns))]

  # --- pseudogeneSize ---
  pseudogene_size <- suppressWarnings(as.integer(cfg[["pseudogeneSize"]]))
  pseudogene_size <- if (length(pseudogene_size) >= 1) pseudogene_size[[1]] else NA_integer_
  if (is.na(pseudogene_size) || pseudogene_size < 1) pseudogene_size <- 2L

  # --- pseudogeneSeed ---
  pseudogene_seed <- suppressWarnings(as.integer(cfg[["pseudogeneSeed"]]))
  pseudogene_seed <- if (length(pseudogene_seed) >= 1) pseudogene_seed[[1]] else NA_integer_
  if (is.na(pseudogene_seed)) pseudogene_seed <- 7L

  # --- pseudogeneControlRegex ---
  pseudogene_control_regex <- cfg[["pseudogeneControlRegex"]]
  if (is.null(pseudogene_control_regex) || length(pseudogene_control_regex) < 1) {
    pseudogene_control_regex <- c("ONE_INTERGENIC_SITE", "NO_SITE")
  }
  pseudogene_control_regex <- as.character(
    unlist(pseudogene_control_regex, use.names = FALSE)
  )
  pseudogene_control_regex <- pseudogene_control_regex[
    !is.na(pseudogene_control_regex) & nzchar(trimws(pseudogene_control_regex))
  ]
  if (length(pseudogene_control_regex) < 1) {
    pseudogene_control_regex <- c("ONE_INTERGENIC_SITE", "NO_SITE")
  }

  # --- lowCountZCutoff ---
  low_count_z_cutoff <- suppressWarnings(as.numeric(cfg[["lowCountZCutoff"]]))
  low_count_z_cutoff <- if (length(low_count_z_cutoff) >= 1) low_count_z_cutoff[[1]] else NA_real_
  if (is.na(low_count_z_cutoff)) low_count_z_cutoff <- -3

  # --- sdCutoff ---
  sd_cutoff <- suppressWarnings(as.numeric(cfg[["sdCutoff"]]))
  sd_cutoff <- if (length(sd_cutoff) >= 1) sd_cutoff[[1]] else NA_real_
  if (is.na(sd_cutoff)) sd_cutoff <- NA_real_

  # --- minGuides ---
  min_guides <- suppressWarnings(as.integer(cfg[["minGuides"]]))
  min_guides <- if (length(min_guides) >= 1) min_guides[[1]] else NA_integer_
  if (is.na(min_guides) || min_guides < 0) min_guides <- 2L

  # --- includeUnexpressed ---
  include_unexpressed <- cfg[["includeUnexpressed"]]
  if (is.null(include_unexpressed)) {
    include_unexpressed <- character(0)
  } else {
    include_unexpressed <- unique(trimws(as.character(unlist(include_unexpressed, use.names = FALSE))))
    include_unexpressed <- include_unexpressed[!is.na(include_unexpressed) & nzchar(include_unexpressed)]
  }

  list(
    guideNegativeControlPatterns = guide_negative_control_patterns,
    guideNegativeControlFuzzy = guide_negative_control_fuzzy,
    pseudogeneSize = pseudogene_size,
    pseudogeneSeed = pseudogene_seed,
    pseudogeneControlRegex = pseudogene_control_regex,
    lowCountZCutoff = low_count_z_cutoff,
    sdCutoff = sd_cutoff,
    minGuides = min_guides,
    includeUnexpressed = include_unexpressed
  )
}


.pu_obs_counts_attach_reference <- function(df, reference_tbl) {
  ref_dedup <- reference_tbl %>%
    distinct(.data$construct_barcode, .data$gene_symbol, .data$gene_id)

  df %>%
    left_join(ref_dedup, by = "construct_barcode")
}

.pu_obs_counts_lognorm <- function(counts_tbl, sample_cols) {
  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] lognorm sample_cols={paste(sample_cols, collapse=', ')}"
  ))
  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] lognorm available_cols={paste(colnames(counts_tbl), collapse=', ')}"
  ))

  missing_cols <- setdiff(sample_cols, colnames(counts_tbl))
  if (length(missing_cols) > 0) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] counts table missing requested analysis sample columns: ",
      "{paste(missing_cols, collapse=', ')}"
    ))
  }

  counts_tbl %>%
    select(any_of(c("construct_barcode", "construct_id")), all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "reads"
    ) %>%
    mutate(
      sample = as.character(.data$sample),
      reads = suppressWarnings(as.numeric(.data$reads))
    ) %>%
    group_by(.data$sample) %>%
    mutate(total_reads = sum(.data$reads, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      lognorm = log2((.data$reads / .data$total_reads) * 1e6 + 1)
    )
}

.pu_obs_counts_remove_filtered_barcodes <- function(lognorms_tbl, ref_sample, z_cutoff = -3) {
  ref_rows <- lognorms_tbl %>%
    filter(.data$sample == ref_sample)

  if (nrow(ref_rows) < 2) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] reference sample {ref_sample} has too few rows after lognorm"
    ))
  }

  ref_sd <- stats::sd(ref_rows$lognorm, na.rm = TRUE)
  if (!is.finite(ref_sd) || ref_sd <= 0) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] reference sample {ref_sample} has non-finite or zero lognorm sd"
    ))
  }

  filtered_barcodes <- ref_rows %>%
    mutate(z_score = as.numeric(scale(.data$lognorm))) %>%
    filter(.data$z_score < z_cutoff) %>%
    pull(.data$construct_barcode) %>%
    unique()

  list(
    filtered_barcodes = filtered_barcodes,
    cleaned = lognorms_tbl %>% filter(!(.data$construct_barcode %in% filtered_barcodes))
  )
}

.pu_obs_counts_calculate_lfcs <- function(cleaned_lognorms, ref_sample) {
  if (!(ref_sample %in% unique(cleaned_lognorms$sample))) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] reference sample {ref_sample} not found in cleaned lognorm table"
    ))
  }

  wide_lognorms <- cleaned_lognorms %>%
    select(.data$construct_barcode, .data$sample, .data$lognorm) %>%
    pivot_wider(names_from = "sample", values_from = "lognorm")

  non_ref_cols <- setdiff(colnames(wide_lognorms), c("construct_barcode", ref_sample))

  wide_lfcs <- wide_lognorms %>%
    mutate(across(all_of(non_ref_cols), ~ .x - .data[[ref_sample]])) %>%
    select(.data$construct_barcode, all_of(non_ref_cols)) %>%
    pivot_longer(
      cols = all_of(non_ref_cols),
      names_to = "sample",
      values_to = "lfc"
    )

  cleaned_lognorms %>%
    left_join(wide_lfcs, by = c("construct_barcode", "sample"))
}

.pu_obs_counts_collapse_replicates <- function(lfcs_tbl, ref_sample = NULL) {

    ref_sample_norm <- if (!is.null(ref_sample) && length(ref_sample) >= 1) {
      x <- as.character(ref_sample[[1]])
      x <- .pu_obs_counts_clean_colname(x)[1]
      if (is.na(x) || !nzchar(x)) NA_character_ else x
    } else {
      NA_character_
    }

  lfcs_tbl %>%
    mutate(
      sample = .pu_obs_counts_clean_colname(as.character(.data$sample)),
      sample_norm = .data$sample,
      is_reference_sample = case_when(
        !is.na(ref_sample_norm) & .data$sample_norm == ref_sample_norm ~ TRUE,
        str_detect(.data$sample_norm, "_t0$") ~ TRUE,
        TRUE ~ FALSE
      ),
      biological_sample = .pu_obs_counts_clean_colname(
        str_replace(.data$sample, "_[Rr][0-9]+$", "")
      )
    ) %>%
    filter(!.data$is_reference_sample) %>%
    group_by(.data$biological_sample, .data$construct_barcode) %>%
    summarize(
      n_replicates = n_distinct(.data$sample),
      avg_lfc = mean(.data$lfc, na.rm = TRUE),
      var_lfc = ifelse(n() > 1, stats::var(.data$lfc, na.rm = TRUE), NA_real_),
      sd_lfc = ifelse(n() > 1, stats::sd(.data$lfc, na.rm = TRUE), NA_real_),
      avg_reads = mean(.data$reads, na.rm = TRUE),
      var_reads = ifelse(n() > 1, stats::var(.data$reads, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    rename(sample = .data$biological_sample)
}


.pu_obs_counts_build_replicate_level_lfcs <- function(lfcs_tbl, ref_sample = NULL) {
  ref_sample_norm <- if (!is.null(ref_sample) && length(ref_sample) >= 1) {
    x <- as.character(ref_sample[[1]])
    x <- .pu_obs_counts_clean_colname(x)[1]
    if (is.na(x) || !nzchar(x)) NA_character_ else x
  } else {
    NA_character_
  }

  lfcs_tbl %>%
    mutate(
      sample = .pu_obs_counts_clean_colname(as.character(.data$sample)),
      sample_norm = .data$sample,
      is_reference_sample = case_when(
        !is.na(ref_sample_norm) & .data$sample_norm == ref_sample_norm ~ TRUE,
        str_detect(.data$sample_norm, "_t0$") ~ TRUE,
        TRUE ~ FALSE
      ),
      biological_sample = .pu_obs_counts_clean_colname(
        str_replace(.data$sample, "_[Rr][0-9]+$", "")
      ),
      replicate_sample = .data$sample
    ) %>%
    filter(!.data$is_reference_sample) %>%
    group_by(.data$biological_sample, .data$replicate_sample, .data$construct_barcode) %>%
    summarize(
      n_replicates = 1L,
      avg_lfc = mean(.data$lfc, na.rm = TRUE),
      var_lfc = NA_real_,
      sd_lfc = NA_real_,
      avg_reads = mean(.data$reads, na.rm = TRUE),
      var_reads = NA_real_,
      .groups = "drop"
    ) %>%
    rename(sample = .data$replicate_sample)
}

.pu_obs_counts_summarize_replicate_positive_probability <- function(
  lfcs_tbl,
  ref_sample,
  reference_tbl,
  new_reference_tbl,
  neg_controls_pattern,
  neg_controls_fuzzy,
  pos_controls,
  neg_controls,
  include_unexpressed = character(0),
  sd_cutoff = NA_real_,
  min_guides = 2L,
  perturbation_tags = "ko_"
) {
  rep_lfcs <- .pu_obs_counts_build_replicate_level_lfcs(
    lfcs_tbl = lfcs_tbl,
    ref_sample = ref_sample
  )

  if (nrow(rep_lfcs) < 1) {
    return(tibble(
      sample = character(0),
      normalizedPerturbation = character(0),
      obs_reps_positive_probability_mean = numeric(0),
      obs_reps_positive_probability_sd = numeric(0),
      obs_reps_positive_probability_n = integer(0),
      obs_reps_positive_probability_status = character(0),
      obs_reps_positive_probability_message = character(0)
    ))
  }

  rep_target_tbl <- .pu_obs_counts_get_targets_table(
    lfcs_collapsed = rep_lfcs,
    reference_tbl = reference_tbl,
    new_reference_tbl = new_reference_tbl,
    neg_controls_pattern = neg_controls_pattern,
    neg_controls_fuzzy = neg_controls_fuzzy,
    pos_controls = pos_controls,
    neg_controls = neg_controls,
    include_unexpressed = include_unexpressed,
    sd_cutoff = sd_cutoff,
    min_guides = min_guides,
    perturbation_tags = perturbation_tags
  ) %>%
    mutate(
      sample = .pu_obs_counts_clean_colname(as.character(.data$sample)),
      biological_sample = .pu_obs_counts_clean_colname(
        str_replace(.data$sample, "_[Rr][0-9]+$", "")
      )
    )

  rep_valid_tbl <- rep_target_tbl %>%
    filter(is.finite(.data$positive_probability))

  rep_summary_tbl <- rep_valid_tbl %>%
    group_by(.data$biological_sample, .data$normalizedPerturbation) %>%
    summarize(
      obs_reps_positive_probability_mean = mean(.data$positive_probability, na.rm = TRUE),
      obs_reps_positive_probability_sd = ifelse(
        dplyr::n() > 1,
        stats::sd(.data$positive_probability, na.rm = TRUE),
        NA_real_
      ),
      obs_reps_positive_probability_n = as.integer(dplyr::n()),
      .groups = "drop"
    ) %>%
    rename(sample = .data$biological_sample)

  rep_status_tbl <- rep_target_tbl %>%
    group_by(.data$biological_sample, .data$normalizedPerturbation) %>%
    summarize(
      n_rows_total = dplyr::n(),
      n_rows_fit_ok = sum(.data$positive_probability_model_status == "fit_ok", na.rm = TRUE),
      n_rows_with_value = sum(is.finite(.data$positive_probability), na.rm = TRUE),
      status_values = paste(
        sort(unique(as.character(.data$positive_probability_model_status[!is.na(.data$positive_probability_model_status)]))),
        collapse = "|"
      ),
      message_values = paste(
        unique(as.character(.data$positive_probability_model_message[
          !is.na(.data$positive_probability_model_message) &
            nzchar(trimws(.data$positive_probability_model_message))
        ])),
        collapse = " | "
      ),
      .groups = "drop"
    ) %>%
    mutate(
      obs_reps_positive_probability_status = dplyr::case_when(
        .data$n_rows_with_value >= 2 ~ "fit_ok_replicates",
        .data$n_rows_with_value == 1 ~ "single_replicate_only",
        .data$n_rows_fit_ok >= 1 ~ "fit_ok_no_numeric_value",
        TRUE ~ "not_fit"
      ),
      obs_reps_positive_probability_message = dplyr::case_when(
        .data$n_rows_with_value >= 2 ~ NA_character_,
        .data$n_rows_with_value == 1 ~ "Only one replicate produced a finite positive_probability.",
        nzchar(.data$message_values) ~ .data$message_values,
        nzchar(.data$status_values) ~ paste("Replicate model statuses:", .data$status_values),
        TRUE ~ "No replicate-specific positive_probability values were produced."
      )
    ) %>%
    transmute(
      sample = .data$biological_sample,
      normalizedPerturbation = .data$normalizedPerturbation,
      obs_reps_positive_probability_status = .data$obs_reps_positive_probability_status,
      obs_reps_positive_probability_message = .data$obs_reps_positive_probability_message
    )

  rep_status_tbl %>%
    left_join(rep_summary_tbl, by = c("sample", "normalizedPerturbation")) %>%
    mutate(
      obs_reps_positive_probability_n = dplyr::coalesce(
        suppressWarnings(as.integer(.data$obs_reps_positive_probability_n)),
        0L
      )
    )
}



.pu_obs_counts_get_neg_controls <- function(reference_tbl, pattern, fuzzy = FALSE) {
  if (isTRUE(fuzzy)) {
    reference_tbl %>%
      filter(str_detect(.data$gene_symbol, pattern)) %>%
      pull(.data$construct_barcode) %>%
      unique()
  } else {
    reference_tbl %>%
      filter(.data$gene_symbol %in% pattern) %>%
      pull(.data$construct_barcode) %>%
      unique()
  }
}

.pu_obs_counts_get_construct_z_scores <- function(
  lfcs_tbl,
  reference_tbl,
  neg_controls_pattern,
  neg_controls_fuzzy
) {
  negs <- purrr::map2(
    neg_controls_pattern,
    neg_controls_fuzzy,
    ~ .pu_obs_counts_get_neg_controls(
      reference_tbl = reference_tbl,
      pattern = .x,
      fuzzy = .y
    )
  ) %>%
    unlist(use.names = FALSE) %>%
    unique()

  neg_df <- lfcs_tbl %>%
    filter(!is.na(.data$avg_lfc)) %>%
    group_by(.data$sample) %>%
    group_split() %>%
    purrr::map(function(x) {
      neg_x <- x %>% filter(.data$construct_barcode %in% negs)
      if (nrow(neg_x) < 2) {
        return(tibble(
          sample = .pu_obs_first_nonempty(x$sample),
          neg_lfc_mean = NA_real_,
          neg_lfc_stdev = NA_real_
        ))
      }
      tibble(
        sample = .pu_obs_first_nonempty(x$sample),
        neg_lfc_mean = mean(neg_x$avg_lfc, na.rm = TRUE),
        neg_lfc_stdev = stats::sd(neg_x$avg_lfc, na.rm = TRUE)
      )
    }) %>%
    bind_rows()


  lfcs_tbl %>%
    left_join(neg_df, by = "sample") %>%
    filter(!is.na(.data$avg_lfc)) %>%
    mutate(
      z = dplyr::case_when(
        !is.finite(.data$neg_lfc_mean) ~ NA_real_,
        !is.finite(.data$neg_lfc_stdev) ~ NA_real_,
        .data$neg_lfc_stdev <= 0 ~ NA_real_,
        TRUE ~ (.data$avg_lfc - .data$neg_lfc_mean) / .data$neg_lfc_stdev
      )
    )
}

.pu_obs_counts_group_pseudogenes <- function(
  annotations,
  pseudogene_size,
  gene_col = "gene_symbol",
  control_regex,
  seed = 7
) {
  remapped_annotations <- annotations
  genes <- remapped_annotations[[gene_col]]
  control_remap <- list()

  for (regex in control_regex) {
    control_genes <- genes[grep(regex, genes)]
    if (length(control_genes) < 1) next

    set.seed(seed)
    control_genes <- sample(control_genes)
    n_controls <- length(control_genes)

    denom <- ceiling(n_controls / pseudogene_size)
    if (denom < 1) denom <- 1L

    for (i in seq_len(n_controls)) {
      gene <- control_genes[[i]]
      gene_number <- i %% denom
      if (gene_number == 0) gene_number <- denom
      control_remap[[gene]] <- paste(regex, "_", as.integer(gene_number), sep = "")
    }
  }

  if (length(control_remap) < 1) {
    return(remapped_annotations)
  }

  final_remap <- tibble(
    gene_symbol = names(control_remap),
    remap = as.character(unlist(control_remap, use.names = FALSE))
  )

  remapped_annotations %>%
    left_join(final_remap, by = "gene_symbol") %>%
    mutate(
      gene_symbol = if_else(is.na(.data$remap), .data$gene_symbol, .data$remap)
    ) %>%
    select(-.data$remap)
}

.pu_obs_counts_calculate_target_lfcs <- function(lfcs_tbl, min_guides = 0L) {
  lfcs_tbl %>%
    group_by(.data$sample, .data$gene_symbol) %>%
    summarize(
      n_guides = n_distinct(.data$construct_barcode),
      mean_target_lfc = mean(.data$avg_lfc, na.rm = TRUE),
      target_z = sum(.data$z, na.rm = TRUE) / sqrt(n_guides),
      target_z_pvalue = 2 * stats::pnorm(abs(.data$target_z), lower.tail = FALSE),
      .groups = "drop"
    ) %>%
    filter(.data$n_guides >= min_guides) %>%
    arrange(.data$target_z_pvalue) %>%
    mutate(
      target_z_fdr = stats::p.adjust(.data$target_z_pvalue, method = "fdr")
    ) %>%
    arrange(desc(.data$n_guides))
}

.pu_obs_counts_get_targets_table <- function(
  lfcs_collapsed,
  reference_tbl,
  new_reference_tbl,
  neg_controls_pattern,
  neg_controls_fuzzy,
  pos_controls,
  neg_controls,
  include_unexpressed = character(0),
  sd_cutoff = NA_real_,
  min_guides = 2L,
  perturbation_tags = "ko_"
) {
  if (is.finite(sd_cutoff)) {
    lfcs_collapsed_filtered <- lfcs_collapsed %>%
      filter(is.na(.data$sd_lfc) | .data$sd_lfc < sd_cutoff)
  } else {
    lfcs_collapsed_filtered <- lfcs_collapsed
  }

  result <- lfcs_collapsed_filtered %>%
    .pu_obs_counts_get_construct_z_scores(
      reference_tbl = new_reference_tbl,
      neg_controls_pattern = neg_controls_pattern,
      neg_controls_fuzzy = neg_controls_fuzzy
    ) %>%
    .pu_obs_counts_attach_reference(reference_tbl = new_reference_tbl) %>%
    .pu_obs_counts_calculate_target_lfcs(min_guides = min_guides) %>%
    arrange(.data$sample) %>%
    mutate(
      perturbation = as.character(.data$gene_symbol),
      normalizedPerturbation = .pu_obs_normalize_perturbation(
        .data$gene_symbol,
        response_set = "crispr",
        perturbation_tags = perturbation_tags
      )
    )

  neg_controls_effective <- unique(c(neg_controls, include_unexpressed))

  result <- result %>%
    .pu_obs_add_control_annotations(
      controls = list(
        positive = .pu_obs_normalize_perturbation(
          pos_controls,
          response_set = "crispr",
          perturbation_tags = perturbation_tags
        ),
        negative = .pu_obs_normalize_perturbation(
          neg_controls_effective,
          response_set = "crispr",
          perturbation_tags = perturbation_tags
        )
      )
    ) %>%
    .pu_obs_add_scaled_target_lfc() %>%
    .pu_obs_fit_positive_probability()
    

  result
}


# -----------------------------
# Exported raw-counts pipeline
# -----------------------------

#' Process a CRISPR observations run from raw guide counts.
#'
#' This is the raw-counts path:
#' raw counts -> lognorm -> low-count filtering on T0 ->
#' guide LFCs -> collapsed replicate LFCs -> guide z-scores ->
#' target summaries -> scaled target LFC -> positive probability.
#'
#' @export
powerup_process_observations_from_counts <- function(
  counts_path,
  reference_path,
  schema_path,
  predictions_path,
  out_observations_dir,
  job_id,
  observation_run_id,
  response_set,
  data_version,
  reference_sample,
  analysis_samples,
  positive_controls_path = NULL,
  negative_controls_path = NULL,
  config_path = NULL,
  seed = 1L,
  perturbation_tags = "ko_"
) {
  if (!identical(tolower(response_set), "crispr")) {
    stop("[powerup][OBSERVATIONS_COUNTS] raw-counts observations mode currently supports response_set='crispr' only")
  }

  powerup_dir_create(out_observations_dir)

  .pu_assert_file_exists(counts_path, "counts_path")
  .pu_assert_file_exists(reference_path, "reference_path")
  .pu_assert_file_exists(schema_path, "schema_path")
  .pu_assert_file_exists(predictions_path, "predictions_path")

  if (!is.null(positive_controls_path) && length(positive_controls_path) >= 1) {
    positive_controls_path <- as.character(positive_controls_path[[1]])
    if (!is.na(positive_controls_path) && nzchar(trimws(positive_controls_path))) {
      .pu_assert_file_exists(positive_controls_path, "positive_controls_path")
    } else {
      positive_controls_path <- NULL
    }
  } else {
    positive_controls_path <- NULL
  }

  if (!is.null(negative_controls_path) && length(negative_controls_path) >= 1) {
    negative_controls_path <- as.character(negative_controls_path[[1]])
    if (!is.na(negative_controls_path) && nzchar(trimws(negative_controls_path))) {
      .pu_assert_file_exists(negative_controls_path, "negative_controls_path")
    } else {
      negative_controls_path <- NULL
    }
  } else {
    negative_controls_path <- NULL
  }

  if (!is.null(config_path) && length(config_path) >= 1) {
    config_path <- as.character(config_path[[1]])
    if (!is.na(config_path) && nzchar(trimws(config_path))) {
      .pu_assert_file_exists(config_path, "config_path")
    } else {
      config_path <- NULL
    }
  } else {
    config_path <- NULL
  }


  reference_sample_raw <- trimws(as.character(reference_sample))[1]
  if (is.na(reference_sample_raw) || !nzchar(reference_sample_raw)) {
    stop("[powerup][OBSERVATIONS_COUNTS] reference_sample must be a non-empty string")
  }

  analysis_samples_raw <- unique(trimws(as.character(analysis_samples)))
  analysis_samples_raw <- analysis_samples_raw[!is.na(analysis_samples_raw) & nzchar(analysis_samples_raw)]

  if (length(analysis_samples_raw) < 2) {
    stop("[powerup][OBSERVATIONS_COUNTS] analysis_samples must include at least the T0 reference and one analysis sample")
  }

  if (!(reference_sample_raw %in% analysis_samples_raw)) {
    stop("[powerup][OBSERVATIONS_COUNTS] reference_sample must be included in analysis_samples")
  }

  reference_sample <- .pu_obs_counts_clean_colname(reference_sample_raw)[1]
  analysis_samples <- unique(.pu_obs_counts_clean_colname(analysis_samples_raw))
  analysis_samples <- analysis_samples[!is.na(analysis_samples) & nzchar(analysis_samples)]

  if (length(analysis_samples) < 2) {
    stop("[powerup][OBSERVATIONS_COUNTS] analysis_samples resolved to fewer than 2 cleaned sample names")
  }

  if (!(reference_sample %in% analysis_samples)) {
    stop(glue(
      "[powerup][OBSERVATIONS_COUNTS] cleaned reference_sample '{reference_sample}' must be included in cleaned analysis_samples: ",
      "{paste(analysis_samples, collapse=', ')}"
    ))
  }

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "reference_sample_raw={reference_sample_raw} reference_sample_clean={reference_sample}"
  ))
  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "analysis_samples_raw={paste(analysis_samples_raw, collapse=', ')}"
  ))
  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "analysis_samples_clean={paste(analysis_samples, collapse=', ')}"
  ))


  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "config_path={if (is.null(config_path)) '<NULL>' else config_path}"
  ))

  cfg <- tryCatch(
    .pu_obs_counts_parse_config(config_path),
    error = function(e) {
      message(glue(
        "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
        "CONFIG_PARSE_FAILED err={conditionMessage(e)}"
      ))
      stop(e)
    }
  )
  schema_obj <- .pu_obs_read_json_file(schema_path)

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "OBSERVATIONS_COUNTS start response_set={response_set} data_version={data_version} seed={seed}"
  ))

  counts_tbl <- .pu_obs_counts_read_counts_required(counts_path)
  message(glue(
  "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
  "counts_tbl columns after normalization={paste(colnames(counts_tbl), collapse=', ')}"
))

  reference_tbl <- .pu_obs_counts_read_reference_required(reference_path)
  message(glue(
  "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
  "reference_tbl columns after normalization={paste(colnames(reference_tbl), collapse=', ')}"
))

  pred_tbl <- .pu_obs_load_predictions_required(predictions_path)

  # Resolve user-supplied control lists exactly, then normalize to CRISPR perturbations later.
  pos_controls_raw <- .pu_obs_read_control_file(positive_controls_path)
  neg_controls_raw <- .pu_obs_read_control_file(negative_controls_path)

  if (length(pos_controls_raw) < 1) {
    stop("[powerup][OBSERVATIONS_COUNTS] positive_controls_path must provide at least one positive control gene")
  }
  if (length(neg_controls_raw) < 1) {
    stop("[powerup][OBSERVATIONS_COUNTS] negative_controls_path must provide at least one negative control gene")
  }

  pos_controls <- unique(trimws(as.character(pos_controls_raw)))
  neg_controls <- unique(trimws(as.character(neg_controls_raw)))

  lognorms_tbl <- .pu_obs_counts_lognorm(
    counts_tbl = counts_tbl,
    sample_cols = analysis_samples
  )

  filter_result <- .pu_obs_counts_remove_filtered_barcodes(
    lognorms_tbl = lognorms_tbl,
    ref_sample = reference_sample,
    z_cutoff = cfg$lowCountZCutoff
  )
  cleaned_lognorms <- filter_result$cleaned
  filtered_barcodes <- filter_result$filtered_barcodes

  lfcs_tbl <- .pu_obs_counts_calculate_lfcs(
    cleaned_lognorms = cleaned_lognorms,
    ref_sample = reference_sample
  )

  lfcs_collapsed <- .pu_obs_counts_collapse_replicates(
    lfcs_tbl = lfcs_tbl,
    ref_sample = reference_sample
  )

  if (nrow(lfcs_collapsed) < 1) {
    stop("[powerup][OBSERVATIONS_COUNTS] no non-reference replicate-collapsed rows remained after collapsing analysis samples")
  }

  if (length(unique(lfcs_collapsed$sample)) < 1) {
    stop("[powerup][OBSERVATIONS_COUNTS] no non-reference biological samples remained after collapsing replicates")
  }

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "collapsed_samples={paste(sort(unique(lfcs_collapsed$sample)), collapse=', ')}"
  ))


  new_reference_tbl <- .pu_obs_counts_step(
    "group_pseudogenes",
    .pu_obs_counts_group_pseudogenes(
      annotations = reference_tbl,
      pseudogene_size = cfg$pseudogeneSize,
      gene_col = "gene_symbol",
      control_regex = cfg$pseudogeneControlRegex,
      seed = cfg$pseudogeneSeed
    )
  )

  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] new_reference_tbl rows={nrow(new_reference_tbl)} cols={paste(colnames(new_reference_tbl), collapse=', ')}"
  ))

  target_tbl <- .pu_obs_counts_step(
    "get_targets_table",
    .pu_obs_counts_get_targets_table(
      lfcs_collapsed = lfcs_collapsed,
      reference_tbl = reference_tbl,
      new_reference_tbl = new_reference_tbl,
      neg_controls_pattern = cfg$guideNegativeControlPatterns,
      neg_controls_fuzzy = cfg$guideNegativeControlFuzzy,
      pos_controls = pos_controls,
      neg_controls = neg_controls,
      include_unexpressed = cfg$includeUnexpressed,
      sd_cutoff = cfg$sdCutoff,
      min_guides = cfg$minGuides,
      perturbation_tags = perturbation_tags
    )
  )


  obs_reps_tbl <- .pu_obs_counts_step(
    "summarize_replicate_positive_probability",
    .pu_obs_counts_summarize_replicate_positive_probability(
      lfcs_tbl = lfcs_tbl,
      ref_sample = reference_sample,
      reference_tbl = reference_tbl,
      new_reference_tbl = new_reference_tbl,
      neg_controls_pattern = cfg$guideNegativeControlPatterns,
      neg_controls_fuzzy = cfg$guideNegativeControlFuzzy,
      pos_controls = pos_controls,
      neg_controls = neg_controls,
      include_unexpressed = cfg$includeUnexpressed,
      sd_cutoff = cfg$sdCutoff,
      min_guides = cfg$minGuides,
      perturbation_tags = perturbation_tags
    )
  )

  target_tbl <- target_tbl %>%
    left_join(
      obs_reps_tbl,
      by = c("sample", "normalizedPerturbation")
    )


  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] target_tbl_pre_mutate rows={nrow(target_tbl)} cols={paste(colnames(target_tbl), collapse=', ')}"
  ))

target_tbl <- .pu_obs_counts_step(
  "decorate_target_tbl",
  target_tbl %>%
    mutate(
      sample = .pu_obs_counts_clean_colname(.data$sample),
      sample_base = .pu_obs_counts_base_sample(.data$sample),
      perturbation = as.character(.data$perturbation),
      normalizedPerturbation = as.character(.data$normalizedPerturbation),
      control = case_when(
        .data$positive_control ~ "positive",
        .data$negative_control ~ "negative",
        TRUE ~ NA_character_
      ),
      category_input = NA_character_,
      observationType = "z_scored_avg_lfc",
      observationValue = suppressWarnings(as.numeric(.data$target_z)),
      primaryObservationType = "z_scored_avg_lfc",
      primaryObservationValue = suppressWarnings(as.numeric(.data$target_z))
    ) %>%
    select(
      .data$sample,
      .data$sample_base,
      .data$perturbation,
      .data$normalizedPerturbation,
      .data$gene_symbol,
      .data$mean_target_lfc,
      .data$target_z,
      .data$target_z_pvalue,
      .data$target_z_fdr,
      .data$n_guides,
      .data$control,
      .data$category_input,
      .data$observationType,
      .data$observationValue,
      .data$primaryObservationType,
      .data$primaryObservationValue,
      .data$control_class,
      .data$positive_control,
      .data$negative_control,
      .data$positive,
      .data$neg_lfc_median,
      .data$pos_lfc_median,
      .data$scaled_target_lfc,
      .data$positive_probability,
      .data$positive_prediction,
      .data$positive_probability_model_status,
      .data$positive_probability_model_message,
      .data$obs_reps_positive_probability_mean,
      .data$obs_reps_positive_probability_sd,
      .data$obs_reps_positive_probability_n,
      .data$obs_reps_positive_probability_status,
      .data$obs_reps_positive_probability_message
    )
)


  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] target_tbl_post_mutate rows={nrow(target_tbl)} cols={paste(colnames(target_tbl), collapse=', ')}"
  ))

  target_tbl <- .pu_obs_counts_step(
    "fit_kde_likelihood",
    .pu_obs_fit_kde_likelihood_by_sample(
      target_tbl = target_tbl,
      obs_col = "target_z",
      bw = "nrd0",
      adjust = 1,
      density_floor = 1e-12
    )
  )

  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] target_tbl_post_kde rows={nrow(target_tbl)} cols={paste(colnames(target_tbl), collapse=', ')}"
  ))

  joined_tbl <- .pu_obs_counts_step(
    "build_joined_tbl",
    .pu_obs_build_joined_tbl(
      target_tbl = target_tbl,
      pred_tbl = pred_tbl,
      response_set = response_set,
      perturbation_tags = perturbation_tags
    )
  )

  message(glue(
    "[powerup][OBSERVATIONS_COUNTS] joined_tbl rows={nrow(joined_tbl)} cols={paste(colnames(joined_tbl), collapse=', ')}"
  ))

  pred_norm_tbl <- pred_tbl %>%
    mutate(
      normalizedPerturbationObsJoin = .pu_obs_normalize_perturbation(
        .data$perturbation,
        response_set = response_set,
        perturbation_tags = perturbation_tags
      )
    )

  observed_perts <- sort(unique(as.character(
    target_tbl$normalizedPerturbation[!is.na(target_tbl$normalizedPerturbation)]
  )))
  predicted_perts <- sort(unique(as.character(
    pred_norm_tbl$normalizedPerturbationObsJoin[!is.na(pred_norm_tbl$normalizedPerturbationObsJoin)]
  )))
  overlapping_perts <- intersect(observed_perts, predicted_perts)

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

  nObsRepsPositiveProbabilityRows <- sum(
    is.finite(suppressWarnings(as.numeric(joined_tbl$obs_reps_positive_probability_mean))),
    na.rm = TRUE
  )

  nObsRepsPositiveProbabilitySdRows <- sum(
    is.finite(suppressWarnings(as.numeric(joined_tbl$obs_reps_positive_probability_sd))),
    na.rm = TRUE
  )


  summary <- list(
    schemaVersion = 4,
    jobId = job_id,
    observationRunId = observation_run_id,
    ok = TRUE,
    mode = "OBSERVATIONS_COUNTS",
    generatedAt = as.character(Sys.time()),
    counts = list(
      nRawCountRows = nrow(counts_tbl),
      nReferenceRows = nrow(reference_tbl),
      nFilteredBarcodes = length(filtered_barcodes),
      nLognormRows = nrow(lognorms_tbl),
      nCleanedLognormRows = nrow(cleaned_lognorms),
      nCollapsedGuideRows = nrow(lfcs_collapsed),
      nTargetLevelRows = nrow(target_tbl),
      nPredictionRows = nrow(pred_tbl),
      nMatchedRows = nrow(matched_rows_only),
      nPriorGaussianRows = nPriorGaussianRows,
      nObsLikelihoodRows = nObsLikelihoodRows,
      nContinuousPosteriorRows = nContinuousPosteriorRows,
      nObsRepsPositiveProbabilityRows = nObsRepsPositiveProbabilityRows,
      nObsRepsPositiveProbabilitySdRows = nObsRepsPositiveProbabilitySdRows,
      nObservedSamples = length(observed_samples),
      nPredictedSamples = length(predicted_samples),
      nOverlappingSamplesExact = length(overlapping_samples),
      nOverlappingSamplesBase = length(overlapping_sample_bases),
      nObservedPerturbations = length(observed_perts),
      nPredictedPerturbations = length(predicted_perts),
      nOverlappingPerturbations = length(overlapping_perts)
    ),
    rawCountsConfig = list(
      referenceSample = reference_sample,
      analysisSamples = analysis_samples,
      guideNegativeControlPatterns = cfg$guideNegativeControlPatterns,
      guideNegativeControlFuzzy = as.list(cfg$guideNegativeControlFuzzy),
      pseudogeneSize = cfg$pseudogeneSize,
      pseudogeneSeed = cfg$pseudogeneSeed,
      pseudogeneControlRegex = cfg$pseudogeneControlRegex,
      lowCountZCutoff = cfg$lowCountZCutoff,
      sdCutoff = if (is.finite(cfg$sdCutoff)) cfg$sdCutoff else NULL,
      minGuides = cfg$minGuides,
      includeUnexpressedCount = length(cfg$includeUnexpressed)
    ),
    controls = list(
      positiveCount = length(pos_controls),
      negativeCount = length(neg_controls),
      includeUnexpressedCount = length(cfg$includeUnexpressed)
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
    schemaVersion = 4,
    jobId = job_id,
    observationRunId = observation_run_id,
    mode = "OBSERVATIONS_COUNTS",
    status = "SUCCEEDED",
    responseSet = response_set,
    dataVersion = data_version,
    generatedAt = as.character(Sys.time()),
    inputs = list(
      countsPath = basename(counts_path),
      referencePath = basename(reference_path),
      schemaPath = basename(schema_path),
      positiveControlsPath = if (!is.null(positive_controls_path) && nzchar(trimws(positive_controls_path))) basename(positive_controls_path) else NULL,
      negativeControlsPath = if (!is.null(negative_controls_path) && nzchar(trimws(negative_controls_path))) basename(negative_controls_path) else NULL,
      predictionsPath = basename(predictions_path),
      configPath = if (!is.null(config_path) && nzchar(trimws(config_path))) basename(config_path) else NULL
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

  # In counts mode, target_tbl already includes positive_probability and related
  # columns, so we intentionally reuse target_tbl for the positiveProbability artifact.
  .pu_obs_write_outputs(
    target_tbl = target_tbl,
    positive_tbl = target_tbl,
    joined_tbl = joined_tbl,
    out_observations_dir = out_observations_dir,
    manifest = manifest,
    summary = summary,
    schema_obj = schema_obj
  )

  message(glue(
    "[powerup][jobId={job_id}][observationRunId={observation_run_id}] ",
    "OBSERVATIONS_COUNTS done target_rows={nrow(target_tbl)} ",
    "prediction_rows={nrow(pred_tbl)} ",
    "matched_rows={nrow(matched_rows_only)} ",
    "filtered_barcodes={length(filtered_barcodes)}"
  ))

  invisible(TRUE)
}