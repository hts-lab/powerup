#' Create a gene expression matrix to predict dependencies
#'
#' This function creates a gene expression matrix from a list of files, each containing bulk gene expression profiles. It returns a samples-by-genes matrix.
#' @param scrna_folder Folder that hosts the Seurat object.
#' @param scrna_name Name of the Seurat object (without .rds).
#' @param scrna_group_cells_by The name of the column in the metadata to group by (default = "seurat_clusters").
#' @param predictions_table The predictions table
#' @keywords seurat predictions
#' @import tidyr
#' @import dplyr
#' @import readr
#' @import stringr
#' @import purrr
#' @import Seurat
#' @export
#' @examples
#' append_predictions_to_seurat("./data","my_seurat_object")
append_predictions_to_seurat <- function(scrna_folder, scrna_name, predictions_table,
                                         scrna_group_cells_by = "seurat_clusters"){
  
  show_msg("Loading {scrna_name} ..")
  scrna_data <- readRDS(glue::glue("{scrna_folder}/{scrna_name}.rds"))
  
  show_msg("Appending predictions ..")
  scrna_data@meta.data <- scrna_data@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(cluster_name = paste0("cluster_",Get_Double_Digit_Number(get(scrna_group_cells_by)))) %>%
    left_join(predictions_table %>% rownames_to_column("cluster_name"), by = "cluster_name") %>%
    column_to_rownames("cell_id")
  
  show_msg("Done")
  return(scrna_data)
  
}