#' Retrieve the training data
#'
#' This function generates the full training data for a given model
#' @param model A model object generated with make_xgb_models
#' @param data The dataset that was passed to make_xgb_models
#' @keywords data
#' @import glue 
#' @export
#' @examples
#' get_original_data(my_models[1],my_data)
get_original_data <- function(model, data){
  
  returned_data <- data[model$sample_names,model$feature_names]
  returned_data <- returned_data %>% mutate(y_value = data[model$sample_names,model$perturbation])
  
  return(returned_data)
  
}