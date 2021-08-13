#' Prepare the lineage features table
#'
#' This function creates a table with one-hot encoded lineage features
#' @param ccle_data A data frame with cell line names as rownames and at least one column matching lineage_column_name.
#' @param our_data A data frame with cell line names as rownames and at least one column matching lineage_column_name.
#' @param lineage_column_name The name of the column that holds the lineage information (default = "lineage").
#' @keywords lineage features
#' @export
#' @examples
#' prep_features_lineage(sample.info)
prep_features_lineage <- function(ccle_data, our_data = NULL, lineage_column_name = "lineage"){
  
  lineage_info <- ccle_data %>%
    dplyr::select(all_of(lineage_column_name))
  
  colnames(lineage_info) <- "lineage"
  
  lineage_info <- lineage_info %>%
    dplyr::mutate(lineage = if_else(lineage == "", "unknown", lineage)) %>%
    rownames_to_column("sample") %>% dplyr::distinct()
  
  if(!is.null(our_data)) lineage_info <- lineage_info %>% dplyr::bind_rows(our_data %>% rownames_to_column("sample") %>% dplyr::distinct())
  
  lineage_info <- lineage_info %>%
    dplyr::mutate(value = 1L) %>% 
    tidyr::pivot_wider("sample", names_from = "lineage", values_from = "value", values_fill = 0L) %>%
    dplyr::rename_with(.cols = -"sample", .fn = ~ paste0("Lin_", .x)) 
  
  
  return(lineage_info)
  
}
