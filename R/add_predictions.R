# add_predictions.R

#' Use a model to make predictions on new data
#'
#' This function takes new data, and produces predictions, error estimates, and SHAP values
#' @param model The model object.
#' @param name The name of the perturbation.
#' @param indx Integer index used, for progress report.
#' @param total Integer of the total number of perturbations passed to this function, for progress report.
#' @param new_data A dataframe of new cases with predictors as columns. Sample names are row names.
#' @keywords model predictions
#' @import xgboost purrr
#' @export
#' @examples
#' make_new_data_predictions(my_model,"ko_ctnnb1",1,1,my_new_data)
make_new_data_predictions <- function(model, name, indx, total, new_data){

  # -----------------------------
  # helpers
  # -----------------------------

  ensure_rownames <- function(df) {
    rn <- rownames(df)
    if (is.null(rn) || length(rn) != nrow(df) || any(is.na(rn)) || any(rn == "")) {
      rownames(df) <- paste0("sample_", seq_len(nrow(df)))
    }
    df
  }

  coerce_numeric_matrix <- function(df) {
    for (nm in colnames(df)) {
      x <- df[[nm]]
      if (is.factor(x)) x <- as.character(x)

      if (is.character(x)) {
        suppressWarnings(x_num <- as.numeric(x))
        if (all(is.na(x_num)) && any(!is.na(x))) {
          stop(glue::glue("new_data column '{nm}' is character and cannot be safely converted to numeric."))
        }
        df[[nm]] <- x_num
      } else if (!is.numeric(x) && !is.integer(x)) {
        suppressWarnings(x_num <- as.numeric(x))
        if (all(is.na(x_num)) && any(!is.na(x))) {
          stop(glue::glue("new_data column '{nm}' has type {class(x)[1]} and cannot be safely converted to numeric."))
        }
        df[[nm]] <- x_num
      }
    }
    m <- as.matrix(df)
    storage.mode(m) <- "double"
    m
  }

  get_xgb_feature_names <- function(booster) {
    attrs <- tryCatch(xgboost::xgb.attributes(booster), error = function(e) NULL)
    if (!is.null(attrs) && "feature_names" %in% names(attrs)) {
      fn <- attrs[["feature_names"]]
      if (!is.null(fn) && length(fn) > 0) return(fn)
    }

    if (!is.null(booster$feature_names) && length(booster$feature_names) > 0) {
      return(booster$feature_names)
    }

    stop("Unable to determine feature names from model$model (xgboost booster).")
  }

  get_xgb_shap_predcontrib <- function(booster, X_mat, sample_names) {
    dm <- xgboost::xgb.DMatrix(data = X_mat)

    phis <- predict(booster, dm, predcontrib = TRUE)
    phis <- as.matrix(phis)

    cn <- colnames(phis)
    bias_idx <- NA_integer_

    if (!is.null(cn)) {
      if ("BIAS" %in% cn) {
        bias_idx <- which(cn == "BIAS")[1]
      } else if ("(Intercept)" %in% cn) {
        bias_idx <- which(cn == "(Intercept)")[1]
      }
    }

    if (is.na(bias_idx) || length(bias_idx) == 0) {
      bias_idx <- ncol(phis)
    }

    bias <- phis[, bias_idx]
    shap <- phis[, -bias_idx, drop = FALSE]

    if (!is.null(colnames(X_mat))) {
      colnames(shap) <- colnames(X_mat)
    }

    shap_df <- as.data.frame(shap, check.names = FALSE)
    rownames(shap_df) <- sample_names

    vals <- apply(shap, 2, function(x) sum(abs(x), na.rm = TRUE))
    contrib <- tibble::tibble(
      term = colnames(shap),
      value = as.numeric(vals)
    ) %>% dplyr::arrange(dplyr::desc(.data$value))

    pos_terms <- contrib %>% dplyr::filter(.data$value > 0) %>% dplyr::pull(.data$term)

    list(
      shap_values = shap_df,
      shap_bias = bias,
      shap_table = contrib,
      good_terms = pos_terms
    )
  }

  # -----------------------------
  # start
  # -----------------------------
  cat(glue::glue("[{lubridate::now('America/New_York')}] Making predictions for {name} ({indx} of {total}) .."), sep = "\n")
  flush.console()

  if (!is.data.frame(new_data)) new_data <- as.data.frame(new_data, check.names = FALSE)
  new_data <- ensure_rownames(new_data)

  if (is.null(model) || is.null(model$model) || !inherits(model$model, "xgb.Booster")) {
    stop("model$model must be an xgboost Booster (xgb.Booster).")
  }
  if (is.null(model$error_model) || !inherits(model$error_model, "xgb.Booster")) {
    stop("model$error_model must be an xgboost Booster (xgb.Booster).")
  }

  feat <- NULL
  if (!is.null(model$feature_names) && length(model$feature_names) > 0) {
    feat <- model$feature_names
  } else {
    feat <- get_xgb_feature_names(model$model)
  }

  missing <- setdiff(feat, colnames(new_data))
  if (length(missing) > 0) {
    stop(glue::glue("new_data is missing {length(missing)} required features, e.g. {missing[[1]]}"))
  }
  new_data <- new_data[, feat, drop = FALSE]

  X_mat <- coerce_numeric_matrix(new_data)
  new_data_dm <- xgboost::xgb.DMatrix(X_mat)

  # -----------------------------
  # Mean prediction
  # -----------------------------
  pred_mean <- as.numeric(stats::predict(model$model, new_data_dm))

  # -----------------------------
  # Variance prediction
  # error_model predicts log( squared_error + eps )
  # -----------------------------
  variance_epsilon <- if (!is.null(model$variance_epsilon)) as.numeric(model$variance_epsilon) else 1e-8
  min_variance <- if (!is.null(model$min_variance)) as.numeric(model$min_variance) else 1e-8
  z95 <- 1.959963984540054

  pred_log_var <- as.numeric(stats::predict(model$error_model, new_data_dm))
  pred_var <- pmax(exp(pred_log_var), min_variance)
  pred_sd <- sqrt(pred_var)

  pred_pi_lower_95 <- pred_mean - z95 * pred_sd
  pred_pi_upper_95 <- pred_mean + z95 * pred_sd

  # -----------------------------
  # Threshold-based event probabilities
  # -----------------------------
  response_cutoff <- if (!is.null(model$response_cutoff)) as.numeric(model$response_cutoff) else NA_real_
  decreasing <- if (!is.null(model$decreasing)) isTRUE(model$decreasing) else FALSE

  prob_below_cutoff <- rep(NA_real_, length(pred_mean))
  prob_above_cutoff <- rep(NA_real_, length(pred_mean))
  prob_target_event <- rep(NA_real_, length(pred_mean))
  target_event_definition <- NA_character_

  if (!is.na(response_cutoff)) {
    standardized <- (response_cutoff - pred_mean) / pred_sd
    prob_below_cutoff <- stats::pnorm(standardized)
    prob_above_cutoff <- 1 - prob_below_cutoff

    if (decreasing) {
      prob_target_event <- prob_below_cutoff
      target_event_definition <- "P(y <= cutoff)"
    } else {
      prob_target_event <- prob_above_cutoff
      target_event_definition <- "P(y >= cutoff)"
    }
  }

  sample_names <- rownames(new_data)
  names(pred_mean) <- sample_names
  names(pred_var) <- sample_names
  names(pred_sd) <- sample_names
  names(pred_log_var) <- sample_names
  names(pred_pi_lower_95) <- sample_names
  names(pred_pi_upper_95) <- sample_names
  names(prob_below_cutoff) <- sample_names
  names(prob_above_cutoff) <- sample_names
  names(prob_target_event) <- sample_names

  # Explain predictions
  shap <- get_xgb_shap_predcontrib(model$model, X_mat, sample_names)

  shap_values_df <- shap$shap_values
  if (!is.data.frame(shap_values_df)) {
    shap_values_df <- as.data.frame(shap_values_df, check.names = FALSE)
  }

  if (is.null(rownames(shap_values_df)) || length(rownames(shap_values_df)) != length(sample_names)) {
    rownames(shap_values_df) <- sample_names
  }

  shap_values_df[["(Intercept)"]] <- as.numeric(shap$shap_bias)

  cols <- colnames(shap_values_df)
  cols <- c("(Intercept)", setdiff(cols, "(Intercept)"))
  shap_values_df <- shap_values_df[, cols, drop = FALSE]


  if (is.null(model$new_data) || !is.list(model$new_data)) {
    model$new_data <- list()
  }
  model$new_data$shap_values <- shap_values_df

  # Attach new data outputs
  model$new_data$data <- new_data

  # legacy compatibility
  model$new_data$predictions <- pred_mean
  model$new_data$predictions_error <- pred_sd

  # new richer outputs
  model$new_data$pred_mean <- pred_mean
  model$new_data$pred_var <- pred_var
  model$new_data$pred_sd <- pred_sd
  model$new_data$pred_log_var <- pred_log_var
  model$new_data$pred_pi_lower_95 <- pred_pi_lower_95
  model$new_data$pred_pi_upper_95 <- pred_pi_upper_95
  model$new_data$prob_below_cutoff <- prob_below_cutoff
  model$new_data$prob_above_cutoff <- prob_above_cutoff
  model$new_data$prob_target_event <- prob_target_event
  model$new_data$target_event_definition <- target_event_definition
  model$new_data$response_cutoff <- response_cutoff
  model$new_data$decreasing <- decreasing

  model$new_data$feature_contribution <- shap$shap_table
  model$new_data$important_features <- shap$good_terms
  model$new_data$shap_bias <- shap$shap_bias

  return(model)
}



#' Use a batch of models to make predictions on new data
#'
#' This function takes a list of models and makes predictions on new data.
#' @param models A list with model objects generated by make_xgb_models.
#' @param new_data A dataframe of new cases with predictors as columns. Sample names are row names.
#' @param models_to_use Optional vector with subset of names of models to use.
#' @keywords model predictions
#' @import xgboost purrr
#' @export
#' @examples
#' add_predictions(models, my_new_data)
add_predictions <- function(models, new_data, models_to_use = NULL){

  # Subset to only needed models if provided
  if (!is.null(models_to_use) && length(models_to_use) > 0) {
    models_to_use <- as.character(models_to_use)
    models_to_use <- models_to_use[nzchar(models_to_use)]

    missing_models <- setdiff(models_to_use, names(models))
    if (length(missing_models) > 0) {
      stop(glue::glue(
        "models_to_use contains {length(missing_models)} model(s) not present in models, e.g. {missing_models[[1]]}"
      ))
    }

    models <- models[models_to_use]
  }


  inputs <- list(
    model = models,
    name  = names(models),
    indx  = seq_along(models)
  )

  models_with_predictions <- purrr::pmap(
    .l = inputs,
    .f = make_new_data_predictions,
    total = length(inputs$indx),
    new_data = new_data
  )

  return(models_with_predictions)
}