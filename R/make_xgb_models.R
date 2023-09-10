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
#' @import xgboost purrr fastshap
#' @export
#' @examples
#' make_xgb_model("ko_ctnnb1",1,1,my_data)
make_xgb_model <- function(perturbation, indx, total, dataset, 
                            response_cutoff = 0.75,
                            weight_cap = 0.05,
                            nfolds = 3, 
                            nrepeats = 3, 
                            nrounds = 100, 
                            min_score = 0.5, 
                            skip_eval = FALSE,
                            shuffle = FALSE,
                            n_threads = 4,
                            use_gpu = TRUE, gpu_id = 0){
  
  cat(glue::glue("[{lubridate::now('US/Eastern')}] Training a model for {perturbation} ({indx} of {total}) .."))
  flush.console()
  
  
  # This keeps one column of dependency scores (renamed 'y_value') plus all predictors
  prepare_model_data <- function(perturbation, data, tag = "ko_", response_cutoff = 0.75, nfolds = 3, nrepeats = 3){
    
    prepared_data <- data %>%
      mutate(y_value = get(perturbation)) %>% 
      select(-starts_with(tag)) %>%
      na.omit() %>% as_tibble(rownames = "cell_line")  
    
    # Create a response column to help stratify cases
    prepared_data <- prepared_data %>% column_to_rownames("cell_line") %>% mutate(response = y_value > response_cutoff)
    
    # Note: We are using the full data as there is no tuning right now.
    data_folds <- vfold_cv(prepared_data, v = nfolds, strata = y_value, repeats = nrepeats, breaks = 20, pool = 0.05)
    
    output <- list()
    output$original_data <- prepared_data
    output$dfolds <- data_folds
    
    return(output)
    
  }
  
  
  # This creates an object that stores model parameters
  # Ideally this is tuning the parameters but we skip this for now
  prepare_model_params <- function(data){
    
    params <- list()
    params$booster <- "gbtree"
    params$objective <- "reg:squarederror"
    
    
    # These parameters seem to do OK on average for all perturbations
    params$eta <- 0.04 # 0.04
    params$gamma <- 0 # 0.01
    params$alpha <- 0.35
    params$lambda <- 0.7
    params$max_depth = 3
    params$subsample = 1
    params$colsample_bytree = 1
    params$colsample_bylevel = 0.2 # 0.2
    params$colsample_bynode = 0.8 # 0.8
    
    
    return(params)
    
  }
  
  # This calculates weights per case based on response
  get_weights <- function(y_value, response_cutoff, weight_cap = 0.05){
    
    # If all samples are above/below the cutoff, return equal weights
    if (weight_cap == 0 || sum(y_value >= response_cutoff) == 0 || sum(y_value < response_cutoff) == 0){
      return(rep(1/length(y_value),times=length(y_value)))
    }
    
    # We count how many cases we have of each response 'status'
    status_counts <- table(if_else(y_value >= response_cutoff, "A", "B"))
    
    # Decide the majority group
    status_major_count <- if_else(status_counts["A"] > status_counts["B"], status_counts["A"], status_counts["B"])
    
    # We calculate the weight of an individual sensitivity status
    status_weight <- 1/status_counts
    
    # We assign the weight to each observation based on observed sensitivity
    weights <- if_else(y_value > response_cutoff, status_weight["A"], status_weight["B"])
    
    # We normalize so that total weights add up to 1
    weights <- weights/sum(weights)
    
    # We cap each observation's individual weight at weight_cap
    weights <- if_else(weights > weight_cap, weight_cap, weights)
    
    # We redistribute the 'lost' weight from the capping step to the remaining samples
    leftover_weight <- (1 - sum(weights))/status_major_count
    weights <- if_else(weights == weight_cap, weight_cap, weights + leftover_weight)
    
    return(weights)
    
  }
  
  # This resamples the given data slice     
  get_weighted_set <- function(data, weights){
    
    data <- data %>%
      slice_sample(prop = 1,
                   replace = TRUE, 
                   weight_by = weights
      )
    
    
    return(data)
    
  }
  
  # This puts the data in DMatrix format for xgboost
  # Optional: To generate a null model we can shuffle the outcome here
  get_DMatrix <- function(data, weights = NULL, shuffle = FALSE){
    
    x_values <- data %>% select(-"y_value",-"response") %>% as.matrix()
    y_values <- data %>% pull(y_value)
    if (shuffle) y_values <- sample(y_values)
    if(!is.null(weights)){
      data <- xgb.DMatrix(data = x_values, label = y_values, weight = 1000*weights) 
    } else {
      data <- xgb.DMatrix(data = x_values, label = y_values) 
    }
    
    
    return(data)
  }
  
  # This calculates SHAP values
  get_xgb_shap <- function(model, data){
    
    pfun <- function(object, newdata) {
      predict(object, newdata = newdata)
    }
    
    shap_obj <- fastshap::explain(model, 
                                  exact = TRUE,
                                  X = data %>% select(-"y_value",-"response") %>% as.matrix(), 
                                  pred_wrapper = pfun, 
                                  adjust = TRUE)
    
    contrib <- tibble(
      term = names(shap_obj),
      value = apply(shap_obj, MARGIN = 2, FUN = function(x) sum(abs(x)))
    ) %>% arrange(desc(value))
    
    nonzero_terms <- contrib %>% filter(value > 0) %>% pull(term)
    
    shap_obj <- shap_obj %>% as.data.frame()
    rownames(shap_obj) <- rownames(data)   
    
    shap_output <- list()
    shap_output$shap_values = shap_obj
    shap_output$shap_table = contrib
    shap_output$good_terms = nonzero_terms
    
    return(shap_output)
    
  }
  
  
  # If the SD is zero, cor() will throw an error
  get_pseudo_cor <- function(x, y){
    
    if (sd(x) == 0 | sd(y) == 0){
      
      x[1] = x[1] + 1e-6
      y[1] = y[1] + 1e-6
      
    } else {
      
      # Do nothing
      
    }
    
    return(cor(x,y))
    
    
  }
  
  
  get_rmse <- function(x, y){
    return( sqrt( mean( (x - y)^2  ) ) )
  }
  
  
  # Step 1: Prepare the data (keep only this perturbation's outcome values and split into folds)
  model_data <- prepare_model_data(perturbation = perturbation, 
                                   data = dataset, 
                                   response_cutoff = response_cutoff, 
                                   nfolds = nfolds, nrepeats = nrepeats)
  
  # Step 2: Define parameters
  model_params <- prepare_model_params(data = model_data)
  
  # Step 3: Assess current parameters with repeated k-fold CV
  # Ideally we are tuning hyperparameters here
  if(!skip_eval){
    
    data_splits = model_data$dfolds$splits
    
    # Grab the analysis parts, for each: resample with bias, and create a training DMatrix
    training_sets <- map(data_splits, analysis)
    
    if(weight_cap > 0){    
      training_weights <- training_sets %>% map(pull, y_value) %>% map(get_weights, response_cutoff = response_cutoff, weight_cap = weight_cap)
      # training_matrices <- training_sets %>% map2(training_weights, get_weighted_set) %>% map(get_DMatrix)
      training_matrices <- training_sets %>% map2(training_weights, get_DMatrix, shuffle = shuffle)
    } else {
      training_matrices <- training_sets %>% map(get_DMatrix, shuffle = shuffle)
    }
    
    # Grab the assessment parts, for each: do as above + pull the y_values for later use
    validation_sets <- map(data_splits, assessment)
    
    if(weight_cap > 0){
      validation_weights <- validation_sets %>% map(pull, y_value) %>% map(get_weights, response_cutoff = response_cutoff, weight_cap = weight_cap)
      #  validation_matrices <- validation_sets %>% map2(validation_weights, get_weighted_set) 
      validation_matrices <- validation_sets %>% map2(validation_weights, get_DMatrix) 
    } else {
      validation_matrices <- validation_sets %>% map(get_DMatrix)
    }
    
    validation_y_values <- map(validation_matrices, getinfo, "label") #instead of pull and "y_value")
    #validation_matrices <- map(validation_matrices, get_DMatrix)
    
    # Use the analysis subsets for creating a model, then use the assessment subset to make predictions and calculate correlation
    score_models <- map(training_matrices, xgboost, 
                        params = model_params,
                        nthread = n_threads,
                        max_bin = 64,
                        tree_method = if_else(use_gpu,"gpu_hist","auto"),
                        gpu_id = gpu_id,
                        nrounds = nrounds,
                        early_stopping_rounds = 10, 
                        verbose = 0)
    
    score_predictions <-  score_models %>% map2(validation_matrices, predict)
    
    scores <- score_predictions %>% map2(validation_y_values, get_pseudo_cor) %>% unlist()
    
    scores_rmse <- score_predictions %>% map2(validation_y_values, get_rmse) %>% unlist()
    
    
  } else {
    
    scores <- rep(1,9)
    scores_rmse <- rep(0,9)
    
  }        
  
  
  cat(glue::glue(" R^2 = {round(mean(scores^2),3)} +/- {round(1.96*sd(scores^2),3)} , RMSE = {round(mean(scores_rmse),5)} , (n={length(scores)})"))
  flush.console()
  
  
  # Prepare output
  output <- list()  
  output$perturbation_name <- perturbation
  output$scores <- scores
  output$scores_rmse <- scores_rmse
  
  # If the score is good enough, we proceed with extra steps                                                                              
  if (!is.na(mean(scores)) & mean(scores^2) >= min_score){
    
    # Fit one last model using all data    
    last_params <- model_params # Ideally we have found the best params and we set them here
    last_nrounds <- nrounds # Ideally this has been tuned too
    
    last_weights <- model_data$original_data %>% pull(y_value) %>% get_weights(response_cutoff = response_cutoff, weight_cap = weight_cap)
    
    if (weight_cap > 0){
      #   last_matrix <-  get_weighted_set(model_data$original_data, last_weights) %>% get_DMatrix()
      last_matrix <-  get_DMatrix(model_data$original_data, last_weights)
    } else {   
      last_matrix <- get_DMatrix(model_data$original_data)
    }
    
    
    # We create a DMatrix using the original (non-resampled) data
    last_validation <- get_DMatrix(model_data$original_data)
    
    # We fit a last model
    last_model <- xgboost(data = last_matrix, 
                          params = last_params, 
                          nthread = n_threads,
                          max_bin = 64,
                          tree_method = if_else(use_gpu,"gpu_hist","auto"),
                          gpu_id = gpu_id,
                          nrounds = last_nrounds, 
                          early_stopping_rounds = 10, verbose = 0)
    
    # We collect the predictions on the same data
    last_predictions <- predict(last_model, newdata = last_validation)
    names(last_predictions) <- rownames(model_data$original_data)
    
    null_prediction <- predict(last_model, newdata = last_matrix) %>% mean()
    
    # Create an error estimate
    errors <- (last_predictions - model_data$original_data$y_value)
    names(errors) <- rownames(model_data$original_data)
    
    # Create a matrix of errors vs features
    error_data <- xgb.DMatrix(data = model_data$original_data %>% select(-"y_value",-"response") %>% as.matrix(),
                              label = errors^2)
    
    # Fit a model on error (using default params)
    error_model <- xgboost(data = error_data, 
                           nrounds = 100, early_stopping_rounds = 10, 
                           max_bin = 64,
                           nthread = n_threads,
                           tree_method = if_else(use_gpu,"gpu_hist","auto"),
                           gpu_id = gpu_id,
                           verbose = 0) 
    
    cat(glue::glue(" E = +/- {round(1.96*sqrt(mean(errors^2,na.rm=T)),3)}"), sep = "\n")
    flush.console()
    
    # Get feature contributions                                                                          
    shap <- get_xgb_shap(last_model, model_data$original_data)
    
    
    # Finish preparing outputs
    output$model <- last_model
    output$error_model <- error_model 
    output$null_prediction <- null_prediction
    output$predictions <- last_predictions
    output$predictions_error <- errors
    output$feature_contribution <- shap$shap_table
    output$important_features <- shap$good_terms
    output$shap_values <- shap$shap_values
    output$sample_names <- rownames(model_data$original_data)
    output$feature_names <- setdiff(colnames(model_data$original_data),c("y_value","response"))
    
    # output$data <- model_data$original_data
    
  } else {
    
    # This model isn't good enough, so we save some time and skip this step.
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
#' @import xgboost purrr fastshap
#' @export
#' @examples
#' fit_depmap_models(my_data, c("ko_ctnnb1","ko_myod1"))
fit_depmap_models <- function(depmap_data, models_to_make,
                              response_cutoff = 0.5,
                              weight_cap = 0,
                              nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.5,
                              skip_eval = FALSE, shuffle = FALSE,
                              n_threads = 4,
                              use_gpu = TRUE, gpu_id = 0){
  
  my_models <- map2(
    models_to_make, seq_along(models_to_make), make_xgb_model,  
    total = length(models_to_make),
    dataset = depmap_data,
    response_cutoff = response_cutoff, weight_cap = weight_cap,
    nfolds = nfolds, nrepeats = nrepeats, nrounds = nrounds, min_score = min_score,
    skip_eval = skip_eval, shuffle = shuffle, n_threads = n_threads, use_gpu = use_gpu, gpu_id = gpu_id)
  
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
#' @import xgboost purrr fastshap tidyverse glue lubridate tidymodels rsample
#' @export
#' @examples
#' fit_models_and_save(my_data, c("ko_ctnnb1","ko_myod1"))
fit_models_and_save <- function(perturbs, chunk_indx, 
                                model_dataset, response_cutoff = 0.5,
                                weight_cap = 0,
                                nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.5,
                                skip_eval = FALSE, shuffle = FALSE,
                                n_threads = 4,
                                use_gpu = TRUE, gpu_id = 0, seed = 123, path = NULL){
  
  library(tidyverse)
  library(glue)
  library(purrr)
  library(lubridate)
  library(tidymodels)
  library(rsample)
  library(xgboost)
  library(fastshap)
  library(mixmap)
  
  if (is.null(path)) path = "."
  
  if (!file.exists(glue::glue("{path}/models_chunk_{chunk_indx}.rds"))){
    
    set.seed(seed)
    
    my_models <- fit_depmap_models(depmap_data = model_dataset, 
                                   models_to_make = perturbs, 
                                   response_cutoff = response_cutoff,
                                   weight_cap = weight_cap,
                                   nfolds = nfolds, nrepeats = nrepeats, nrounds = nrounds, min_score = min_score,
                                   skip_eval = skip_eval, shuffle = shuffle,
                                   n_threads = n_threads,
                                   use_gpu = use_gpu, gpu_id = gpu_id)
    
    
    
    saveRDS(my_models,glue::glue("{path}/models_chunk_{chunk_indx}.rds"))
    
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
#' @import xgboost purrr furrr future fastshap tidyverse glue lubridate tidymodels rsample
#' @export
#' @examples
#' fit_models_in_parallel(my_data, c("ko_ctnnb1","ko_myod1"))
fit_models_in_parallel <- function(perturbs, chunk_size = 20, 
                                   model_dataset, response_cutoff = 0.5,
                                   weight_cap = 0,
                                   nfolds = 3, nrepeats = 1, nrounds = 200, min_score = 0.5,
                                   skip_eval = FALSE, shuffle = FALSE, 
                                   n_threads = 4,
                                   use_gpu = TRUE, gpu_id = c(0), seed = 123, path = NULL){

  perturb_splits <- split(perturbs, ceiling(seq_along(perturbs)/chunk_size))
  
  # Generate a list of inputs
  inputs <- list()
  inputs$perturbs <- perturb_splits
  inputs$chunk_indx <- seq_along(perturb_splits)
  inputs$gpu_id <- rep(gpu_id,length.out=length(perturb_splits))
  
  furrr::future_pmap(inputs,fit_models_and_save,
                     model_dataset = model_dataset, response_cutoff = response_cutoff,
                     weight_cap = weight_cap,
                     nfolds = nfolds, nrepeats = nrepeats, nrounds = nrounds, min_score = min_score,
                     skip_eval = skip_eval, shuffle = shuffle,
                     n_threads = n_threads,
                     use_gpu = use_gpu, seed = seed, path = path)
  
  return("Done")
  
}
