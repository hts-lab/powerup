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
    # Prevent as.matrix(data.frame) from turning everything into character
    for (nm in colnames(df)) {
      x <- df[[nm]]
      if (is.factor(x)) x <- as.character(x)

      if (is.character(x)) {
        suppressWarnings(x_num <- as.numeric(x))
        # if there are any non-NA original values but numeric is all NA => unsafe conversion
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

  # Robustly extract feature names from the fitted booster
  get_xgb_feature_names <- function(booster) {
    # Preferred: attribute "feature_names" (most consistent)
    attrs <- tryCatch(xgboost::xgb.attributes(booster), error = function(e) NULL)
    if (!is.null(attrs) && "feature_names" %in% names(attrs)) {
      fn <- attrs[["feature_names"]]
      if (!is.null(fn) && length(fn) > 0) return(fn)
    }

    # Sometimes stored directly
    if (!is.null(booster$feature_names) && length(booster$feature_names) > 0) {
      return(booster$feature_names)
    }

    stop("Unable to determine feature names from model$model (xgboost booster).")
  }

  # ----
  # SHAP for xgboost boosters via predcontrib (exact for trees)
  # ----
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

    # Fallback: last column is bias
    if (is.na(bias_idx) || length(bias_idx) == 0) {
      bias_idx <- ncol(phis)
    }

    bias <- phis[, bias_idx]
    shap <- phis[, -bias_idx, drop = FALSE]

    # Ensure feature names match X
    if (!is.null(colnames(X_mat))) {
      colnames(shap) <- colnames(X_mat)
    }

    shap_df <- as.data.frame(shap, check.names = FALSE)
    rownames(shap_df) <- sample_names

    # feature importance table: sum abs shap per feature
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
  cat(glue::glue("[{lubridate::now('US/Eastern')}] Making predictions for {name} ({indx} of {total}) .."), sep = "\n")
  flush.console()

  if (!is.data.frame(new_data)) new_data <- as.data.frame(new_data, check.names = FALSE)
  new_data <- ensure_rownames(new_data)

  # Keep only the features needed by the model
  if (is.null(model) || is.null(model$model) || !inherits(model$model, "xgb.Booster")) {
    stop("model$model must be an xgboost Booster (xgb.Booster).")
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

  # Convert to numeric matrix + DMatrix
  X_mat <- coerce_numeric_matrix(new_data)
  new_data_dm <- xgboost::xgb.DMatrix(X_mat)

  # Make predictions and error estimates for each sample
  predictions <- as.numeric(stats::predict(model$model, new_data_dm))
  error <- as.numeric(stats::predict(model$error_model, new_data_dm))

  sample_names <- rownames(new_data)
  names(predictions) <- sample_names
  names(error) <- sample_names

  # Explain predictions (BYPASS fastshap completely)
  shap <- get_xgb_shap_predcontrib(model$model, X_mat, sample_names)

  # ------------------------------------------------------------
  # Include SHAP bias as a special "feature" column in shap_values
  # so downstream plotting can reconstruct:
  #   pred ≈ (Intercept) + sum(feature SHAP)
  #
  # We name it exactly "(Intercept)" for consistency with common
  # SHAP conventions. Keep check.names=FALSE so parentheses remain.
  # ------------------------------------------------------------
  shap_values_df <- shap$shap_values
  if (!is.data.frame(shap_values_df)) {
    shap_values_df <- as.data.frame(shap_values_df, check.names = FALSE)
  }

  # Ensure rownames match sample_names (defensive)
  if (is.null(rownames(shap_values_df)) || length(rownames(shap_values_df)) != length(sample_names)) {
    rownames(shap_values_df) <- sample_names
  }

  # Add bias column as "(Intercept)"
  shap_values_df[["(Intercept)"]] <- as.numeric(shap$shap_bias)

  # Deterministic column order: put intercept first, then features in their existing order
  cols <- colnames(shap_values_df)
  cols <- c("(Intercept)", setdiff(cols, "(Intercept)"))
  shap_values_df <- shap_values_df[, cols, drop = FALSE]

  # Replace shap_values stored on model with the augmented version
  model$new_data$shap_values <- shap_values_df



  # Attach new data outputs to the original model
  model$new_data$data <- new_data
  model$new_data$predictions <- predictions
  model$new_data$predictions_error <- error

  model$new_data$feature_contribution <- shap$shap_table
  model$new_data$important_features <- shap$good_terms
  #model$new_data$shap_values <- shap$shap_values (see above)
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
  if(!is.null(models_to_use) && length(models_to_use) > 0) {
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