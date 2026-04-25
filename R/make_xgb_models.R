# make_xgb_models.R


#' Make predictive models of dependencies
#'
#' This function creates an XGBoost model
#' @param perturbation Column name of the perturbation (e.g. "ko_ctnnb1").
#' @param indx Integer index used, for progress report.
#' @param total Integer of the total number of perturbations passed to this function, for progress report.
#' @param dataset A dataframe with the perturbation in a column and all other predictors. Sample names are row names.
#' @param response_cutoff The value above which the sample is considered sensitive.
#' @param weight_cap The maximum weight of each minority case when resampling. Set to 0 if no resampling needed.
#' @param nfolds The number of folds in k-fold cross validation.
#' @param nrepeats The number of repeats in k-fold cross validation.
#' @param nrounds The maximum number of trees in the XGBoost model.
#' @param min_score The minimum number of r value for a model to be considered for the next stage (making predictions and calculating SHAP values).
#' @param skip_eval Default = FALSE. If TRUE, k-fold CV will not be conducted and instead all models will be pushed to the next stage.
#' @param use_gpu Default = TRUE. Set to FALSE if using CPU.
#' @keywords model
#' @import xgboost purrr
#' @import dplyr tibble glue lubridate rsample magrittr
#' @export
#' @examples
#' make_xgb_model("ko_ctnnb1",1,1,my_data)
make_xgb_model <- function(perturbation, indx, total, dataset,
                            response_cutoff = 0.75, decreasing = F,
                            weight_cap = 0.05,
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
                            cor_data = NULL, cor_n_features = 1000,
                            perturbation_tags = "ko_",
                            use_gpu = TRUE, gpu_id = 0){

  cat(glue::glue("[{lubridate::now('America/New_York')}] Training a model for {perturbation} ({indx} of {total}) .."))
  flush.console()

  shap_top_n <- as.integer(Sys.getenv("POWERUP_SHAP_TOP_N", unset = "100"))
  if (is.na(shap_top_n) || shap_top_n < 1) {
    shap_top_n <- 100L
  }

  # Small epsilon for variance modeling
  variance_epsilon <- 1e-8
  min_variance <- 1e-8
  normal_z_95 <- 1.959963984540054

  # Final full-data refit should use a fixed number of rounds chosen from CV,
  # not early stopping against the full training data.
  last_nrounds <- as.integer(nrounds)
  eval_count <- max(1L, as.integer(nfolds) * as.integer(nrepeats))

  # This keeps one column of dependency scores (renamed 'y_value') plus all predictors
  prepare_model_data <- function(
    perturbation,
    data,
    perturbation_tags = "ko_",
    response_cutoff = 0.75,
    nfolds = 3,
    nrepeats = 3,
    cor_data = NULL,
    cor_num = 1000
  ) {

    if (!(perturbation %in% colnames(data))) {
      stop(glue::glue("Perturbation column '{perturbation}' not found in dataset"))
    }

    tags <- as.character(perturbation_tags %||% "ko_")
    tags <- trimws(tags)
    tags <- unique(tags[nzchar(tags)])

    if (length(tags) < 1) {
      tags <- "ko_"
    }

    escape_regex <- function(x) {
      gsub("([][{}()+*^$|\\\\.?-])", "\\\\\\1", x)
    }

    prefix_pattern <- paste0("^(", paste(vapply(tags, escape_regex, character(1)), collapse = "|"), ")")

    drop_prefixed_cols <- function(df) {
      out <- df %>% dplyr::select(-dplyr::any_of(perturbation))
      prefixed_cols <- grep(prefix_pattern, colnames(out), value = TRUE, perl = TRUE)
      if (length(prefixed_cols) > 0) {
        out <- out %>% dplyr::select(-dplyr::any_of(prefixed_cols))
      }
      out
    }

    if (is.null(cor_data)) {

      prepared_data <- data %>%
        dplyr::mutate(y_value = .data[[perturbation]]) %>%
        drop_prefixed_cols() %>%
        stats::na.omit() %>%
        tibble::as_tibble(rownames = "cell_line")

    } else {

      correlated_features <- NULL

      if (perturbation %in% colnames(cor_data)) {
        correlated_features <- cor_data %>%
          dplyr::top_n(cor_num, abs(.data[[perturbation]])) %>%
          dplyr::pull(feature)
      }

      if (!is.null(correlated_features) && length(correlated_features) > 0) {

        prepared_data <- data %>%
          dplyr::mutate(y_value = .data[[perturbation]]) %>%
          drop_prefixed_cols() %>%
          dplyr::select(y_value, dplyr::any_of(correlated_features)) %>%
          stats::na.omit() %>%
          tibble::as_tibble(rownames = "cell_line")

      } else {

        prepared_data <- data %>%
          dplyr::mutate(y_value = .data[[perturbation]]) %>%
          drop_prefixed_cols() %>%
          stats::na.omit() %>%
          tibble::as_tibble(rownames = "cell_line")

      }
    }

    prepared_data <- prepared_data %>%
      tibble::column_to_rownames("cell_line") %>%
      dplyr::mutate(response = y_value >= response_cutoff)

    data_folds <- rsample::vfold_cv(
      prepared_data,
      v = nfolds,
      strata = y_value,
      repeats = nrepeats,
      breaks = 20,
      pool = 0.05
    )

    output <- list()
    output$original_data <- prepared_data
    output$dfolds <- data_folds

    return(output)
  }


  # This creates an object that stores model parameters
  prepare_model_params <- function(data, xgb_params){

    if (is.null(xgb_params)) xgb_params <- list()

    params <- list()
    params$booster <- "gbtree"
    params$objective <- "reg:squarederror"

    params$eta <- 0.04
    params$gamma <- 0
    params$alpha <- 0.35
    params$lambda <- 0.7
    params$sampling_method = "gradient_based"
    params$colsample_bytree = 1
    params$colsample_bylevel = 0.2
    params$colsample_bynode = 0.8

    if (is.null(xgb_params$eta)) params$eta <- 0.04 else params$eta <- xgb_params$eta
    if (is.null(xgb_params$gamma)) params$gamma <- 0 else params$gamma <- xgb_params$gamma
    if (is.null(xgb_params$alpha)) params$alpha <- 0.35 else params$alpha <- xgb_params$alpha
    if (is.null(xgb_params$lambda)) params$lambda <- 0.7 else params$lambda <- xgb_params$lambda
    if (is.null(xgb_params$colsample_bytree)) params$colsample_bytree <- 1 else params$colsample_bytree <- xgb_params$colsample_bytree
    if (is.null(xgb_params$colsample_bylevel)) params$colsample_bylevel <- 0.2 else params$colsample_bylevel <- xgb_params$colsample_bylevel
    if (is.null(xgb_params$colsample_bynode)) params$colsample_bynode <- 0.8 else params$colsample_bynode <- xgb_params$colsample_bynode

    return(params)
  }

  # This calculates weights per case based on response
  get_weights <- function(y_value, response_cutoff, weight_cap = 0.05){

    if (weight_cap == 0 || sum(y_value >= response_cutoff) == 0 || sum(y_value < response_cutoff) == 0){
      return(rep(1/length(y_value), times = length(y_value)))
    }

    status_counts <- table(dplyr::if_else(y_value >= response_cutoff, "A", "B"))
    status_major_count <- dplyr::if_else(status_counts["A"] > status_counts["B"], status_counts["A"], status_counts["B"])
    status_weight <- 1 / status_counts

    weights <- dplyr::if_else(y_value >= response_cutoff, status_weight["A"], status_weight["B"])
    weights <- weights / sum(weights)
    weights <- dplyr::if_else(weights > weight_cap, weight_cap, weights)

    leftover_weight <- (1 - sum(weights)) / status_major_count
    weights <- dplyr::if_else(weights == weight_cap, weight_cap, weights + leftover_weight)

    return(weights)
  }

  # This puts the data in DMatrix format for xgboost
  get_DMatrix <- function(data, weights = NULL, shuffle = FALSE){

    x_df <- data %>% dplyr::select(-"y_value", -"response")
    x_values <- as.matrix(x_df)
    storage.mode(x_values) <- "double"

    y_values <- data %>% dplyr::pull(y_value)
    if (shuffle) y_values <- sample(y_values)

    if (!is.null(weights)) {
      data <- xgboost::xgb.DMatrix(data = x_values, label = y_values, weight = 1000 * weights)
    } else {
      data <- xgboost::xgb.DMatrix(data = x_values, label = y_values)
    }

    return(data)
  }

  # This calculates SHAP values
  get_xgb_shap <- function(model, X, sample_names = NULL) {
    cat("[SHAP] ENTER get_xgb_shap\n")
    flush.console()

    X_mat <- as.matrix(X)
    storage.mode(X_mat) <- "double"

    if (!inherits(model, "xgb.Booster")) {
      stop("get_xgb_shap currently supports xgb.Booster only (fastshap bypass).")
    }

    cat("[SHAP] Using xgboost predcontrib=TRUE (bypass fastshap)\n")
    flush.console()

    dm <- xgboost::xgb.DMatrix(data = X_mat)
    phis <- predict(model, dm, predcontrib = TRUE)
    phis <- as.matrix(phis)

    cn <- colnames(phis)
    bias_idx <- NA_integer_
    if (!is.null(cn)) {
      if ("BIAS" %in% cn) bias_idx <- which(cn == "BIAS")[1]
      if (is.na(bias_idx) && "(Intercept)" %in% cn) bias_idx <- which(cn == "(Intercept)")[1]
    }
    if (is.na(bias_idx)) bias_idx <- ncol(phis)

    bias <- phis[, bias_idx]
    shap <- phis[, -bias_idx, drop = FALSE]

    if (is.null(sample_names)) {
      sample_names <- rownames(X)
    }
    if (!is.null(sample_names)) {
      rownames(shap) <- sample_names
      names(bias) <- sample_names
    }

    if (!is.null(colnames(X_mat))) {
      colnames(shap) <- colnames(X_mat)
    }

    shap_with_bias <- cbind(bias, shap)
    colnames(shap_with_bias)[1] <- "(Intercept)"

    shap_values_df <- data.frame(shap_with_bias, check.names = FALSE)

    if (is.null(sample_names)) sample_names <- rownames(X)
    if (!is.null(sample_names) && length(sample_names) == nrow(shap_values_df)) {
      rownames(shap_values_df) <- sample_names
    }

    cols <- colnames(shap_values_df)
    cols <- c("(Intercept)", setdiff(cols, "(Intercept)"))
    shap_values_df <- shap_values_df[, cols, drop = FALSE]

    vals <- apply(shap, 2, function(x) sum(abs(x), na.rm = TRUE))
    contrib <- tibble::tibble(
      term = colnames(shap),
      value = as.numeric(vals)
    ) %>% dplyr::arrange(dplyr::desc(.data$value))

    pos_terms <- contrib %>% dplyr::filter(.data$value > 0) %>% dplyr::pull(.data$term)

    cat(sprintf("[SHAP] predcontrib phis dim: %d x %d\n", nrow(phis), ncol(phis)))
    cat(sprintf("[SHAP] shap dim (no bias): %d x %d\n", nrow(shap), ncol(shap)))
    cat(sprintf("[SHAP] shap dim (WITH bias): %d x %d\n", nrow(shap_with_bias), ncol(shap_with_bias)))
    cat("[SHAP] EXIT get_xgb_shap OK\n")
    flush.console()

    return(list(
      shap_values = shap_values_df,
      bias = bias,
      shap_table = contrib,
      good_terms = pos_terms
    ))
  }

  .pu_collapse_shap_rowwise_matrix <- function(shap_df, top_n) {

    # Convert rownames → cell_line
    if (!("cell_line" %in% colnames(shap_df))) {
      shap_df <- shap_df %>% tibble::rownames_to_column("cell_line")
    }
    
    feature_cols <- setdiff(colnames(shap_df), "cell_line")

    if (any(feature_cols %in% c("__other__", "__other_count__"))) {
      stop("Reserved SHAP names detected in upstream features")
    }

    has_intercept <- "(Intercept)" %in% feature_cols

    out_rows <- vector("list", nrow(shap_df))

    for (i in seq_len(nrow(shap_df))) {

      row <- shap_df[i, , drop = FALSE]
      cl <- row$cell_line

      feats <- row[, feature_cols, drop = FALSE]

      intercept_val <- NULL
      if (has_intercept) {
        intercept_val <- feats[["(Intercept)"]]
        feats[["(Intercept)"]] <- NULL
      }

      vals <- suppressWarnings(as.numeric(feats[1, ]))
      names(vals) <- colnames(feats)

      vals[is.na(vals)] <- 0

      if (length(vals) <= top_n) {
        kept <- vals
        other_sum <- 0
        other_count <- 0
      } else {
        ord <- order(-abs(vals), names(vals))
        keep_idx <- ord[seq_len(top_n)]

        kept <- vals[keep_idx]

        drop_idx <- setdiff(seq_along(vals), keep_idx)

        other_sum <- sum(vals[drop_idx])
        other_count <- length(drop_idx)
      }

      row_out <- c(
        list(cell_line = cl),
        if (has_intercept) list(`(Intercept)` = intercept_val) else NULL,
        as.list(kept),
        list(
          "__other__" = other_sum,
          "__other_count__" = as.integer(other_count)
        )
      )

      out_rows[[i]] <- row_out
    }

    out_df <- dplyr::bind_rows(out_rows)

    base_cols <- c("cell_line")
    if (has_intercept) base_cols <- c(base_cols, "(Intercept)")

    feature_cols_new <- setdiff(colnames(out_df), c(base_cols, "__other__", "__other_count__"))
    feature_cols_new <- sort(feature_cols_new)

    out_df <- out_df[, c(base_cols, feature_cols_new, "__other__", "__other_count__"), drop = FALSE]

    return(out_df)
  }



  # If SD is zero/NA, correlation is not well-defined
  get_pseudo_cor <- function(x, y){

    x <- as.numeric(x)
    y <- as.numeric(y)

    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]

    if (length(x) < 2 || length(y) < 2) {
      return(NA_real_)
    }

    sx <- stats::sd(x)
    sy <- stats::sd(y)

    if (is.na(sx) || is.na(sy) || sx == 0 || sy == 0) {
      return(NA_real_)
    }

    return(stats::cor(x, y))
  }

  get_rmse <- function(x, y){

    x <- as.numeric(x)
    y <- as.numeric(y)

    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) {
      return(NA_real_)
    }

    return(sqrt(mean((x[ok] - y[ok])^2)))
  }

  get_R2 <- function(x, y){

    x <- as.numeric(x)
    y <- as.numeric(y)

    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]

    if (length(x) < 2 || length(y) < 2) {
      return(NA_real_)
    }

    denom <- sum((y - mean(y))^2)

    if (is.na(denom) || denom <= 0) {
      return(NA_real_)
    }

    return(1 - sum((x - y)^2) / denom)
  }

  get_discrete_sensitivity <- function(pred, obs, discrete_cut, decreasing = F){

    if (decreasing){
      pred_d = dplyr::if_else(pred <= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs <= discrete_cut, T, F)
    } else {
      pred_d = dplyr::if_else(pred >= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs >= discrete_cut, T, F)
    }

    TP = sum(pred_d & obs_d)
    FN = sum(!pred_d & obs_d)

    denom <- TP + FN
    if (denom == 0) return(NA_real_)

    result = TP / denom
    return(result)
  }


  get_discrete_specificity <- function(pred, obs, discrete_cut, decreasing = F){

    if (decreasing){
      pred_d = dplyr::if_else(pred <= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs <= discrete_cut, T, F)
    } else {
      pred_d = dplyr::if_else(pred >= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs >= discrete_cut, T, F)
    }

    TN = sum(!pred_d & !obs_d)
    FP = sum(pred_d & !obs_d)

    denom <- TN + FP
    if (denom == 0) return(NA_real_)

    result = TN / denom
    return(result)
  }


  get_discrete_fpr <- function(pred, obs, discrete_cut, decreasing = F){

    if (decreasing){
      pred_d = dplyr::if_else(pred <= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs <= discrete_cut, T, F)
    } else {
      pred_d = dplyr::if_else(pred >= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs >= discrete_cut, T, F)
    }

    FP = sum(pred_d & !obs_d)
    TN = sum(!pred_d & !obs_d)

    denom <- FP + TN
    if (denom == 0) return(NA_real_)

    result = FP / denom
    return(result)
  }


  get_discrete_ppv <- function(pred, obs, discrete_cut, decreasing = F){

    if (decreasing){
      pred_d = dplyr::if_else(pred <= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs <= discrete_cut, T, F)
    } else {
      pred_d = dplyr::if_else(pred >= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs >= discrete_cut, T, F)
    }

    TP = sum(pred_d & obs_d)
    FP = sum(pred_d & !obs_d)

    denom <- TP + FP
    if (denom == 0) return(NA_real_)

    result = TP / denom
    return(result)
  }


  get_discrete_npv <- function(pred, obs, discrete_cut, decreasing = F){

    if (decreasing){
      pred_d = dplyr::if_else(pred <= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs <= discrete_cut, T, F)
    } else {
      pred_d = dplyr::if_else(pred >= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs >= discrete_cut, T, F)
    }

    TN = sum(!pred_d & !obs_d)
    FN = sum(!pred_d & obs_d)

    denom <- TN + FN
    if (denom == 0) return(NA_real_)

    result = TN / denom
    return(result)
  }


  get_discrete_accuracy <- function(pred, obs, discrete_cut, decreasing = F){

    if (decreasing){
      pred_d = dplyr::if_else(pred <= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs <= discrete_cut, T, F)
    } else {
      pred_d = dplyr::if_else(pred >= discrete_cut, T, F)
      obs_d = dplyr::if_else(obs >= discrete_cut, T, F)
    }

    TP = sum(pred_d & obs_d)
    FP = sum(pred_d & !obs_d)
    TN = sum(!pred_d & !obs_d)
    FN = sum(!pred_d & obs_d)

    denom <- TP + FP + TN + FN
    if (denom == 0) return(NA_real_)

    result = (TP + TN) / denom
    return(result)
  }
  

  # Step 1: Prepare the data
  model_data <- prepare_model_data(
    perturbation = perturbation,
    data = dataset,
    perturbation_tags = perturbation_tags,
    response_cutoff = response_cutoff,
    nfolds = nfolds,
    nrepeats = nrepeats,
    cor_data = cor_data,
    cor_num = cor_n_features
  )

  n_features_used <- ncol(model_data$original_data) - 2L

  # Step 2: Define parameters
  model_params <- prepare_model_params(data = model_data, xgb_params = xgb_params)

  oof_predictions_df <- NULL

  # Step 3: Assess current parameters with repeated k-fold CV
  if (!skip_eval){

    data_splits = model_data$dfolds$splits

    training_sets <- purrr::map(data_splits, rsample::analysis)

    if (weight_cap > 0){
      training_weights <- training_sets %>%
        purrr::map(dplyr::pull, y_value) %>%
        purrr::map(get_weights, response_cutoff = response_cutoff, weight_cap = weight_cap)
      training_matrices <- training_sets %>%
        purrr::map2(training_weights, get_DMatrix, shuffle = shuffle)
    } else {
      training_matrices <- training_sets %>% purrr::map(get_DMatrix, shuffle = shuffle)
    }

    validation_sets <- purrr::map(data_splits, rsample::assessment)

    if (weight_cap > 0){
      validation_weights <- validation_sets %>%
        purrr::map(dplyr::pull, y_value) %>%
        purrr::map(get_weights, response_cutoff = response_cutoff, weight_cap = weight_cap)
      validation_matrices <- validation_sets %>% purrr::map2(validation_weights, get_DMatrix)
    } else {
      validation_matrices <- validation_sets %>% purrr::map(get_DMatrix)
    }

    validation_y_values <- purrr::map(validation_matrices, xgboost::getinfo, "label")
    validation_sample_names <- purrr::map(validation_sets, rownames)

    cv_params <- model_params
    cv_params$max_depth <- max_depth
    cv_params$subsample <- f_subsample
    cv_params$nthread <- n_threads
    cv_params$max_bin <- 64
    cv_params$tree_method <- if (use_gpu) "gpu_hist" else "auto"
    if (use_gpu) cv_params$gpu_id <- gpu_id

    score_models <- purrr::map2(
      training_matrices, validation_matrices,
      function(dtrain, dval) {
        xgboost::xgb.train(
          params = cv_params,
          data = dtrain,
          nrounds = nrounds,
          evals = list(train = dtrain, eval = dval),
          early_stopping_rounds = 10,
          verbose = 0
        )
      }
    )

    cv_best_iterations <- purrr::map_int(
      score_models,
      function(m) {
        bi <- m$best_iteration
        if (is.null(bi) || is.na(bi) || bi < 1) {
          return(as.integer(nrounds))
        }
        as.integer(bi)
      }
    )

    last_nrounds <- max(
      1L,
      as.integer(round(stats::median(cv_best_iterations, na.rm = TRUE)))
    )

    pred_best <- function(model, dmat) {
      bi <- model$best_iteration
      if (!is.null(bi) && !is.na(bi)) {
        return(predict(model, dmat, iterationrange = c(0, bi)))
      }
      predict(model, dmat)
    }

    score_predictions <- purrr::map2(score_models, validation_matrices, pred_best)

    # Build OOF prediction table for variance modeling
    oof_predictions_df <- purrr::pmap_dfr(
      list(
        pred = score_predictions,
        obs = validation_y_values,
        cell_line = validation_sample_names,
        fold_index = seq_along(score_predictions)
      ),
      function(pred, obs, cell_line, fold_index) {
        tibble::tibble(
          cell_line = as.character(cell_line),
          pred = as.numeric(pred),
          obs = as.numeric(obs),
          fold_index = as.integer(fold_index)
        )
      }
    )

    scores <- score_predictions %>% purrr::map2(validation_y_values, get_pseudo_cor) %>% unlist()
    scores_rmse <- score_predictions %>% purrr::map2(validation_y_values, get_rmse) %>% unlist()
    scores_R2 <- score_predictions %>% purrr::map2(validation_y_values, get_R2) %>% unlist()

    scores_d_sensitivity <- score_predictions %>% purrr::map2(validation_y_values, get_discrete_sensitivity, response_cutoff, decreasing) %>% unlist()
    scores_d_specificity <- score_predictions %>% purrr::map2(validation_y_values, get_discrete_specificity, response_cutoff, decreasing) %>% unlist()
    scores_d_fpr <- score_predictions %>% purrr::map2(validation_y_values, get_discrete_fpr, response_cutoff, decreasing) %>% unlist()
    scores_d_ppv <- score_predictions %>% purrr::map2(validation_y_values, get_discrete_ppv, response_cutoff, decreasing) %>% unlist()
    scores_d_npv <- score_predictions %>% purrr::map2(validation_y_values, get_discrete_npv, response_cutoff, decreasing) %>% unlist()
    scores_d_accuracy <- score_predictions %>% purrr::map2(validation_y_values, get_discrete_accuracy, response_cutoff, decreasing) %>% unlist()

    rm(score_models)
    rm(data_splits)
    rm(training_matrices)
    rm(validation_matrices)

  } else {

    scores <- rep(1, eval_count)
    scores_R2 <- rep(1, eval_count)
    scores_rmse <- rep(0, eval_count)
    scores_d_sensitivity <- rep(1, eval_count)
    scores_d_specificity <- rep(1, eval_count)
    scores_d_fpr <- rep(0, eval_count)
    scores_d_ppv <- rep(1, eval_count)
    scores_d_npv <- rep(1, eval_count)
    scores_d_accuracy <- rep(1, eval_count)

    # Fallback pseudo-OOF table if evaluation is skipped
    oof_predictions_df <- tibble::tibble(
      cell_line = rownames(model_data$original_data),
      pred = as.numeric(model_data$original_data$y_value),
      obs = as.numeric(model_data$original_data$y_value),
      fold_index = 1L
    )
  }

  cat(glue::glue(
    " r = {round(mean(scores, na.rm = TRUE),3)} +/- {round(ifelse(sum(is.finite(scores)) > 1, 1.96*stats::sd(scores, na.rm = TRUE)/sqrt(sum(is.finite(scores))), NA_real_),3)}",
    " | R2 = {round(mean(scores_R2, na.rm = TRUE),3)} +/- {round(ifelse(sum(is.finite(scores_R2)) > 1, 1.96*stats::sd(scores_R2, na.rm = TRUE)/sqrt(sum(is.finite(scores_R2))), NA_real_),3)}",
    " | RMSE = {round(mean(scores_rmse, na.rm = TRUE),5)}",
    " | CV best_nrounds median = {last_nrounds}",
    " | (n={length(scores)})"
  ))

  flush.console()

  output <- list()
  output$perturbation_name <- perturbation
  output$n_features_used <- n_features_used
  output$scores <- scores
  output$scores_R2 <- scores_R2
  output$scores_rmse <- scores_rmse
  output$scores_d_sensitivity <- scores_d_sensitivity
  output$scores_d_specificity <- scores_d_specificity
  output$scores_d_fpr <- scores_d_fpr
  output$scores_d_ppv <- scores_d_ppv
  output$scores_d_npv <- scores_d_npv
  output$scores_d_accuracy <- scores_d_accuracy
  output$response_cutoff <- response_cutoff
  output$decreasing <- decreasing
  output$variance_epsilon <- variance_epsilon
  output$min_variance <- min_variance
  output$distribution_family <- "gaussian_oof_log_variance"

  # If the score is good enough, proceed
  mean_score <- mean(scores, na.rm = TRUE)
  mean_r2 <- mean(scores_R2, na.rm = TRUE)

  if (is.finite(mean_score) && is.finite(mean_r2) && mean_r2 >= min_score){
    last_params <- model_params

    last_weights <- model_data$original_data %>%
      dplyr::pull(y_value) %>%
      get_weights(response_cutoff = response_cutoff, weight_cap = weight_cap)

    if (weight_cap > 0){
      last_matrix <- get_DMatrix(model_data$original_data, last_weights, shuffle = shuffle)
    } else {
      last_matrix <- get_DMatrix(model_data$original_data, shuffle = shuffle)
    }

    final_params <- last_params
    final_params$max_depth <- max_depth
    final_params$subsample <- f_subsample
    final_params$nthread <- n_threads
    final_params$max_bin <- 64
    final_params$tree_method <- if (use_gpu) "gpu_hist" else "auto"
    if (use_gpu) final_params$gpu_id <- gpu_id

    last_model <- xgboost::xgb.train(
      params = final_params,
      data = last_matrix,
      nrounds = last_nrounds,
      verbose = 0
    )

    last_predictions <- predict(last_model, newdata = last_matrix)
    names(last_predictions) <- rownames(model_data$original_data)

    # We calculcate null prediction from SHAP bias terms  
    # null_prediction <- predict(last_model, newdata = last_matrix) %>% mean()

    # In-sample residuals retained for debugging only
    errors <- (last_predictions - model_data$original_data$y_value)
    names(errors) <- rownames(model_data$original_data)

    # -----------------------------
    # Variance model from OOF residuals
    # -----------------------------
    if (is.null(oof_predictions_df) || nrow(oof_predictions_df) < 2) {
      stop(glue::glue("OOF predictions were not available to fit uncertainty model for {perturbation}"))
    }

    oof_var_df <- oof_predictions_df %>%
      dplyr::mutate(
        sq_error = (.data$pred - .data$obs)^2,
        log_sq_error = log(.data$sq_error + variance_epsilon)
      ) %>%
      dplyr::group_by(.data$cell_line) %>%
      dplyr::summarise(
        y_value = mean(.data$log_sq_error, na.rm = TRUE),
        response = FALSE,
        .groups = "drop"
      ) %>%
      tibble::column_to_rownames("cell_line")

    var_feature_df <- model_data$original_data %>%
      dplyr::select(-"y_value", -"response") %>%
      tibble::rownames_to_column("cell_line")

    oof_var_df <- oof_var_df %>%
      tibble::rownames_to_column("cell_line") %>%
      dplyr::inner_join(var_feature_df, by = "cell_line") %>%
      tibble::column_to_rownames("cell_line")

    if (nrow(oof_var_df) < 10) {
      stop(glue::glue("Too few OOF rows after joining features for uncertainty model: {nrow(oof_var_df)}"))
    }

    error_data <- get_DMatrix(oof_var_df, weights = NULL, shuffle = FALSE)

    err_params <- last_params
    err_params$max_depth <- max_depth
    err_params$subsample <- f_subsample
    err_params$nthread <- n_threads
    err_params$max_bin <- 64
    err_params$tree_method <- if (use_gpu) "gpu_hist" else "auto"
    if (use_gpu) err_params$gpu_id <- gpu_id

    error_model <- xgboost::xgb.train(
      params = err_params,
      data = error_data,
      nrounds = last_nrounds,
      evals = list(train = error_data),
      verbose = 0
    )

    train_log_var_pred <- as.numeric(predict(error_model, newdata = last_matrix))
    train_pred_var <- pmax(exp(train_log_var_pred), min_variance)
    train_pred_sd <- sqrt(train_pred_var)

    cat(glue::glue(" mean predicted 95% half-width = {round(1.96*mean(train_pred_sd, na.rm = TRUE), 3)}"), sep = "\n")
    flush.console()

    X_feat <- model_data$original_data %>%
      dplyr::select(-"y_value", -"response") %>%
      as.matrix()

    shap <- get_xgb_shap(last_model, X_feat, sample_names = rownames(model_data$original_data))

    null_prediction <- unname(shap$bias[[1]])

    output$model <- last_model
    output$error_model <- error_model
    output$null_prediction <- null_prediction

    # legacy + new training outputs
    output$predictions <- last_predictions
    output$predictions_error <- train_pred_sd
    output$predictions_var <- train_pred_var
    output$predictions_log_var <- train_log_var_pred
    output$predictions_residual <- errors
    output$predictions_pi_lower_95 <- as.numeric(last_predictions - normal_z_95 * train_pred_sd)
    output$predictions_pi_upper_95 <- as.numeric(last_predictions + normal_z_95 * train_pred_sd)

    output$feature_contribution <- shap$shap_table
    output$important_features <- shap$good_terms

    collapsed_shap <- .pu_collapse_shap_rowwise_matrix(
      shap$shap_values,
      top_n = shap_top_n
    )

    output$shap_values <- collapsed_shap
    output$shap_bias <- shap$bias

    output$sample_names <- rownames(model_data$original_data)
    output$feature_names <- colnames(X_feat)

    rm(last_model)
    rm(error_model)
    rm(last_matrix)
    rm(shap)
    gc()

  } else {

    cat(glue::glue(" Skipped"), sep = "\n")
    flush.console()
  }

  return(output)
}


#' Make a list of predictive models of dependencies
#'
#' This function creates an XGBoost model for each perturbation given, and returns a list of model objects.
#' @param perturbation Column name of the perturbation (e.g. "ko_ctnnb1").
#' @param indx Integer index used, for progress report.
#' @param total Integer of the total number of perturbations passed to this function, for progress report.
#' @param dataset A dataframe with the perturbation in a column and all other predictors. Sample names are row names.
#' @param response_cutoff The value above which the sample is considered sensitive.
#' @param weight_cap The maximum weight of each minority case when resampling. Set to 0 if no resampling needed.
#' @param nfolds The number of folds in k-fold cross validation.
#' @param nrepeats The number of repeats in k-fold cross validation.
#' @param nrounds The maximum number of trees in the XGBoost model.
#' @param min_score The minimum number of r^2 value for a model to be considered for the next stage (making predictions and calculating SHAP values).
#' @param skip_eval Default = FALSE. If TRUE, k-fold CV will not be conducted and instead all models will be pushed to the next stage.
#' @param use_gpu Default = TRUE. Set to FALSE if using CPU.
#' @keywords model
#' @import xgboost purrr
#' @import dplyr tibble glue lubridate rsample magrittr
#' @export
#' @examples
#' fit_depmap_models(my_data, c("ko_ctnnb1","ko_myod1"))
fit_depmap_models <- function(depmap_data, models_to_make,
                              perturbation_tags = "ko_",
                              response_cutoff = 0.5, decreasing = FALSE,
                              weight_cap = 0,
                              nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.01,
                              max_depth = 3,
                              f_subsample = 1,
                              skip_eval = FALSE, shuffle = FALSE,
                              n_threads = 4,
                              xgb_params = NULL,
                              cor_data = NULL, cor_n_features = 1000,
                              use_gpu = TRUE, gpu_id = 0){
  
  my_models <- purrr::map2(
    models_to_make, seq_along(models_to_make), make_xgb_model,  
    total = length(models_to_make),
    dataset = depmap_data,
    perturbation_tags = perturbation_tags,
    response_cutoff = response_cutoff, decreasing = decreasing,
    weight_cap = weight_cap,
    nfolds = nfolds, nrepeats = nrepeats, nrounds = nrounds, min_score = min_score,
    max_depth = max_depth,
    f_subsample = f_subsample,
    skip_eval = skip_eval, shuffle = shuffle, 
    xgb_params = xgb_params,
    n_threads = n_threads, 
    cor_data = cor_data, cor_n_features = cor_n_features,
    use_gpu = use_gpu, gpu_id = gpu_id)
  
  names(my_models) <- models_to_make
  
  return(my_models)
  
}




#' Make a list of models and save as file
#'
#' This function creates an XGBoost model for each perturbation given, saves the list of models, and returns a message.
#' @param perturbs A vector of perturbations.
#' @param chunk_indx Integer index used, for progress report.
#' @param model_dataset A dataframe with the perturbation in a column and all other predictors. Sample names are row names.
#' @param response_cutoff The value above which the sample is considered sensitive.
#' @param weight_cap The maximum weight of each minority case when resampling. Set to 0 if no resampling needed.
#' @param nfolds The number of folds in k-fold cross validation.
#' @param nrepeats The number of repeats in k-fold cross validation.
#' @param nrounds The maximum number of trees in the XGBoost model.
#' @param min_score The minimum number of r^2 value for a model to be considered for the next stage (making predictions and calculating SHAP values).
#' @param skip_eval Default = FALSE. If TRUE, k-fold CV will not be conducted and instead all models will be pushed to the next stage.
#' @param use_gpu Default = TRUE. Set to FALSE if using CPU.
#' @param seed Random seed
#' @param path Folder path (e.g. "/home/test/models") to save models in.
#' @keywords model
#' @import xgboost purrr glue lubridate rsample
#' @export
#' @examples
#' fit_models_and_save(my_data, c("ko_ctnnb1","ko_myod1"))
fit_models_and_save <- function(perturbs, chunk_indx, 
                                model_dataset, 
                                perturbation_tags = "ko_",
                                response_cutoff = 0.5, decreasing = FALSE,
                                weight_cap = 0,
                                nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.01,
                                max_depth = 3,
                                f_subsample = 1,
                                skip_eval = FALSE, shuffle = FALSE,
                                xgb_params = NULL,
                                n_threads = 4,
                                cor_data = NULL, cor_n_features = 1000,
                                use_gpu = TRUE, gpu_id = 0, seed = 123, path = NULL){
  
  library(glue)
  library(purrr)
  library(lubridate)
  library(rsample)
  library(xgboost)
  
  if (is.null(path)) path = "."
  
  if (!file.exists(glue::glue("{path}/models_chunk_{chunk_indx}.rds"))){
    
    set.seed(seed)
    
    my_models <- fit_depmap_models(depmap_data = model_dataset, 
                                   models_to_make = perturbs, 
                                   perturbation_tags = perturbation_tags,
                                   response_cutoff = response_cutoff, decreasing = decreasing,
                                   weight_cap = weight_cap,
                                   nfolds = nfolds, nrepeats = nrepeats, nrounds = nrounds,
                                   min_score = min_score,
                                   max_depth = max_depth,
                                   f_subsample = f_subsample,
                                   skip_eval = skip_eval, shuffle = shuffle,
                                   xgb_params = xgb_params,
                                   cor_data = cor_data,
                                   cor_n_features = cor_n_features,
                                   n_threads = n_threads,
                                   use_gpu = use_gpu, gpu_id = gpu_id)
    
    
    
    saveRDS(my_models,glue::glue("{path}/models_chunk_{chunk_indx}.rds"))
    
    # Clean up
    rm(my_models)
    gc()
    
    return(glue::glue("Done chunk {chunk_indx}"))
    
  } else {
    
    
    return(glue::glue("Chunk already done {chunk_indx}"))
    
    
  }
  
  
}




#' Make a list of models and save as file (parallel)
#'
#' This function creates an XGBoost model for each perturbation given, saves the list of models, and returns a message.
#' @param perturbs A vector of perturbations.
#' @param chunk_indx Integer index used, for progress report.
#' @param model_dataset A dataframe with the perturbation in a column and all other predictors. Sample names are row names.
#' @param response_cutoff The value above which the sample is considered sensitive.
#' @param weight_cap The maximum weight of each minority case when resampling. Set to 0 if no resampling needed.
#' @param nfolds The number of folds in k-fold cross validation.
#' @param nrepeats The number of repeats in k-fold cross validation.
#' @param nrounds The maximum number of trees in the XGBoost model.
#' @param min_score The minimum number of r^2 value for a model to be considered for the next stage (making predictions and calculating SHAP values).
#' @param skip_eval Default = FALSE. If TRUE, k-fold CV will not be conducted and instead all models will be pushed to the next stage.
#' @param use_gpu Default = TRUE. Set to FALSE if using CPU.
#' @param seed Random seed
#' @param path Folder path (e.g. "/home/test/models") to save models in.
#' @keywords model
#' @import xgboost future glue lubridate rsample
#' @export
#' @examples
#' fit_models_in_parallel(my_data, c("ko_ctnnb1","ko_myod1"))
fit_models_in_parallel <- function(perturbs, chunk_size = 20, 
                                   model_dataset, 
                                   perturbation_tags = "ko_",
                                   response_cutoff = 0.5, decreasing = FALSE,
                                   weight_cap = 0,
                                   nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.01,
                                   max_depth = 3,
                                   f_subsample = 1,
                                   skip_eval = FALSE, shuffle = FALSE, 
                                   xgb_params = NULL,
                                   cor_data = NULL, cor_n_features = 1000,
                                   n_threads = 4,
                                   use_gpu = TRUE, gpu_id = c(0), seed = 123, path = NULL){

  perturb_splits <- split(perturbs, ceiling(seq_along(perturbs)/chunk_size))
  
  # Generate a list of inputs
  inputs <- list()
  inputs$perturbs <- perturb_splits
  inputs$chunk_indx <- seq_along(perturb_splits)
  inputs$gpu_id <- rep(gpu_id,length.out=length(perturb_splits))
  
  furrr::future_pmap(inputs,fit_models_and_save,
                     model_dataset = model_dataset, 
                     perturbation_tags = perturbation_tags,                     
                     response_cutoff = response_cutoff, decreasing = decreasing,
                     weight_cap = weight_cap,
                     nfolds = nfolds, nrepeats = nrepeats, nrounds = nrounds, 
                     min_score = min_score,
                     max_depth = max_depth,
                     f_subsample = f_subsample,
                     skip_eval = skip_eval, shuffle = shuffle,
                     xgb_params = xgb_params,
                     cor_data = cor_data, cor_n_features = cor_n_features,
                     n_threads = n_threads,
                     use_gpu = use_gpu, seed = seed, path = path,  
                     .options = furrr_options(seed = TRUE))
  
  return("Done")
  
}






#' Make a list of models and save as file (parallel)
#'
#' This function creates an XGBoost model for each perturbation given, saves the list of models, and returns a message.
#' @param perturbs A vector of perturbations.
#' @param chunk_indx Integer index used, for progress report.
#' @param model_dataset A dataframe with the perturbation in a column and all other predictors. Sample names are row names.
#' @param response_cutoff The value above which the sample is considered sensitive.
#' @param weight_cap The maximum weight of each minority case when resampling. Set to 0 if no resampling needed.
#' @param nfolds The number of folds in k-fold cross validation.
#' @param nrepeats The number of repeats in k-fold cross validation.
#' @param nrounds The maximum number of trees in the XGBoost model.
#' @param min_score The minimum number of r^2 value for a model to be considered for the next stage (making predictions and calculating SHAP values).
#' @param skip_eval Default = FALSE. If TRUE, k-fold CV will not be conducted and instead all models will be pushed to the next stage.
#' @param use_gpu Default = TRUE. Set to FALSE if using CPU.
#' @param seed Random seed
#' @param path Folder path (e.g. "/home/test/models") to save models in.
#' @keywords model
#' @import xgboost purrr furrr future glue lubridate rsample
#' @export
#' @examples
#' fit_models(my_data, c("ko_ctnnb1","ko_myod1"))
fit_models <- function(perturbs, model_dataset, splits = 10,
                       chunk_size = 6, 
                                   perturbation_tags = "ko_",   
                                   response_cutoff = 0.5, decreasing = FALSE,
                                   weight_cap = 0,
                                   nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.01,
                                   max_depth = 3,
                                   f_subsample = 1,
                                   skip_eval = FALSE, shuffle = FALSE, 
                                   xgb_params = NULL,
                                   cor_data = NULL, cor_n_features = 1000,
                                   n_threads = 4,
                                   use_gpu = TRUE, gpu_id = c(0), seed = 123, path = NULL){
  
  show_msg(glue::glue("[{lubridate::now('America/New_York')}] We will start fitting models for {length(perturbs)} perturbations.
    Data will be stored at  {path}/results"))
  
  
  big_chunks <- split(perturbs, 1:splits)

  # Create folders to host chunk outputs
  big_chunk_count = 0
  for (this_big_chunk in big_chunks){
    big_chunk_count = big_chunk_count + 1
    dir.create(
      file.path(path, "results", glue::glue("big_chunk_{big_chunk_count}")),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  # Release memory
  future::plan(sequential)
  
  
  # Loop through chunks to train models in parallel
  big_chunk_count = 0
  for (this_big_chunk in big_chunks){
    big_chunk_count = big_chunk_count + 1
    show_msg("[{lubridate::now('America/New_York')}] Processing {big_chunk_count} of {splits} - STARTED")
    future::plan(multisession, workers = 8)
    fit_models_in_parallel(perturbs = this_big_chunk, 
                           perturbation_tags = perturbation_tags,
                           chunk_size = chunk_size,
                           min_score = min_score,
                           max_depth = max_depth,
                           model_dataset = model_dataset, 
                           response_cutoff = response_cutoff,
                           path = glue("{path}/results/big_chunk_{big_chunk_count}"),
                           shuffle = shuffle,
                           n_threads = n_threads,
                           use_gpu = use_gpu, gpu_id = gpu_id)
    future::plan(sequential)
    show_msg("[{lubridate::now('America/New_York')}] Processing {big_chunk_count} of {splits} - DONE")
    
  }
  
  
  
  show_msg(glue::glue("[{lubridate::now('America/New_York')}] Done fitting models in parallel. Merging outputs...
    Data will be stored at  {path}/results"))
  
  
  chunk_files <- list.files(glue("{path}/results"), pattern = "models_chunk_*",full.names = T,recursive =  T)
  
  all_chunks <- list()
  
  for(chunk_file in chunk_files){
    
    this_chunk <- readRDS(chunk_file)
    
    all_chunks <- append(all_chunks, this_chunk)
    
  }
  
  
  show_msg(glue::glue("[{lubridate::now('America/New_York')}] Done fitting models in parallel. Saving...
    Data will be stored at {path}/results/models.rds"))
  
  saveRDS(all_chunks,glue::glue("{path}/results/models.rds"))
  
  show_msg(glue("[{lubridate::now('America/New_York')}] Done. The overall average model accuracy (r) was {all_chunks %>% purrr::map('scores') %>% purrr::map(~mean(.x, na.rm = TRUE)) %>% unlist() %>% mean(na.rm = TRUE)}"))

  return("Done")
  
}