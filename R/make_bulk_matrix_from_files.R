#' Create a gene expression matrix to predict dependencies
#'
#' This function creates a gene expression matrix from a list of files, each containing bulk gene expression profiles. It returns a samples-by-genes matrix.
#' @param bulk_files A vector of filenames for gene counts or pre-calculated TPMs.
#' @param pseudocount An arbitrary small number to add before log-transformation (default = 1).
#' @keywords message
#' @import tidyr
#' @import dplyr
#' @import readr
#' @import stringr
#' @import purrr
#' @import edgeR
#' @export
#' @examples
#' show_msg("The number of cells is {cell_count}")
make_bulk_matrix_from_files <- function(bulk_files, pseudocount = 1){
    
  load_bulk_data <- function(file, prior_count = 1){

    data <- NULL

    if (str_detect(file,"gene_reads.gct")){

        data <- readr::read_tsv(file, skip = 2, col_types = cols())
        data <- data %>% select(-"Name") %>% dplyr::rename("gene" = "Description")

        data <- data %>% dplyr::group_by(gene) %>% dplyr::summarize_all(max)

        data <- data %>% filter(!is.na(gene)) %>% column_to_rownames("gene") %>%
        as.matrix() %>% edgeR::cpm(log = T, prior.count = prior_count)

        data[data < 0] <- 0

        data <- data %>% as_tibble(rownames = "gene")

    } else {

        data <- readr::read_tsv(file, col_types = cols())
        colnames(data) <- stringr::word(colnames(data),1,sep=".bam")

    }

    return(data)

   }  
    
    
    # Load bulk files 
    bulk_data <- purrr::map(rev(bulk_files), load_bulk_data, prior_count = pseudocount) %>% 
    purrr::reduce(left_join, by = "gene") %>% 
    na.omit() %>%
    dplyr::rowwise() %>% dplyr::mutate(var = var(dplyr::c_across(-"gene")))

    # Arrange by variance
    bulk_data <- bulk_data %>% dplyr::filter(var > 0) %>% dplyr::arrange(desc(var)) %>% dplyr::ungroup()

    # Convert to matrix
    bulk_data <- bulk_data %>% dplyr::select(-"var") %>% column_to_rownames("gene") %>% t %>% as.matrix()
    
    return(bulk_data)
    
}