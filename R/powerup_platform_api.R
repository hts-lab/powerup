# powerup_platform_api.R

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
})

# ---- small utilities ----

powerup_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

powerup_write_json <- function(path, obj) {
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE), con = path)
}

# ---- CONTRACT 1: preprocess ----
# This is a thin deterministic builder to wire to existing split logic later.
powerup_preprocess <- function(
  gene_expression_path,
  response_path,
  matrix_path,
  out_preprocess_dir,
  data_version,
  response_set,
  seed,
  job_id
) {
  powerup_dir_create(out_preprocess_dir)

  set.seed(as.integer(seed))

  # TODO: replace these reads/transforms with existing logic
  # gene_expression <- readr::read_csv(gene_expression_path, show_col_types = FALSE)
  # response_df     <- readr::read_csv(response_path, show_col_types = FALSE)
  # user_matrix     <- readr::read_tsv(matrix_path, show_col_types = FALSE)

  # TODO: create deterministic perturbations table with modelKey -> perturbation as source of truth.
  # Must include columns: modelKey + perturbation
  #
  # Write:
  # - perturbations.csv
  # - train_set.csv
  # - test_set.csv
  # - user_samples.csv
  # - manifest.json (must contain seed)

  # Placeholder to force wiring (remove once implemented):
  stop("powerup_preprocess(): wire this to existing preprocessing/splitting logic.")
}

# ---- CONTRACT 2: train models ----
# This function receives a list of model keys and trains them.
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
  weight_cap = 0.05,
  nfolds = 3,
  nrepeats = 3,
  nrounds = 100,
  max_depth = 3,
  f_subsample = 1,
  min_score = 0.5,
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

  # Load tables (or split artifacts) that have already been preprocessed.
  train_df <- readr::read_csv(train_set_path, show_col_types = FALSE) %>% tibble::column_to_rownames("cell_line")
  test_df  <- readr::read_csv(test_set_path,  show_col_types = FALSE) %>% tibble::column_to_rownames("cell_line")
  user_df  <- readr::read_csv(user_samples_path, show_col_types = FALSE) %>% tibble::column_to_rownames("cell_line")

  # perturbations table is the mapping source of truth
  pert_tbl <- readr::read_csv(perturbations_path, show_col_types = FALSE)
  if (!("modelKey" %in% names(pert_tbl))) stop("perturbations.csv must contain column: modelKey")

  # Build quick lookup from modelKey -> perturbation column name
  # EXPECTATION: pert_tbl has a column "perturbation" that matches a column in train_df/test_df
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

  # We train each model and write artifacts per modelKey
  total <- length(model_keys)

  for (i in seq_along(model_keys)) {
    mk <- model_keys[[i]]
    perturbation <- key_map$perturbation[key_map$modelKey == mk][[1]]

    model_out_dir <- file.path(out_models_dir, mk)
    powerup_dir_create(model_out_dir)

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

    # Clean memory between models
    rm(fit)
    gc()
  }

  invisible(TRUE)
}

# ---- CONTRACT 3: finalize ----
powerup_finalize <- function(shard_manifests_dir, out_aggregates_dir, job_id) {
  powerup_dir_create(out_aggregates_dir)

  # TODO: Aggregation logic:
  # - read shard manifests
  # - compile model_status.csv and a final manifest.json
  stop("powerup_finalize(): implement aggregation of shard outputs.")
}