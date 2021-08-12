#' Create a gene expression matrix to predict dependencies
#'
#' This function creates a gene expression matrix from a list of files, each containing bulk gene expression profiles. It returns a samples-by-genes matrix.
#' @param scrna_folder Folder that hosts the Seurat object.
#' @param scrna_name Name of the Seurat object (without .rds).
#' @param scrna_group_cells_by The name of the column in the metadata to group by (default = "seurat_clusters").
#' @param scrna_assay_markers The assay (e.g. "RNA", "SCT", "alra") used for differential expression analysis.
#' @param scrna_assay_vars The assay (e.g. "RNA", "SCT", "alra") used for pre-calculated variable features.
#' @param scrna_test_markers The method used (e.g. "MAST", "wilcox", "t") used for differential expression analysis. See help for FindAllMarkers() in Seurat.
#' @param scrna_logfc See help for FindAllMarkers() in Seurat.
#' @param scrna_minpct See help for FindAllMarkers() in Seurat.
#' @param scrna_label_prefix When creating a pseudobulk profile for each subgroup, this is a prefix to add to their labels (default = "cluster").
#' @param pseudocount An arbitrary small number to add to gene counts before log transformation. 1 works well for 10X data, 0.083 works well for Drop-seq data (default = 1).
#' @param qvalue_cutoff FDR cutoff to consider a marker significant.
#' @param force_find_markers By default, finding markers will be skipped if done previously (saved as [scrna_name].[scrna_group_cells_by].markers inside the [scrna_folder] folder). If TRUE, finding markers will always happen.
#' @keywords pseudobulk
#' @import tidyr
#' @import dplyr
#' @import readr
#' @import stringr
#' @import purrr
#' @import edgeR
#' @import Seurat
#' @export
#' @examples
#' show_msg("The number of cells is {cell_count}")
make_bulk_matrix_from_seurat <- function(scrna_folder,
                                               scrna_name,
                                               scrna_group_cells_by = "seurat_clusters",
                                               scrna_assay_markers = "SCT",
                                               scrna_assay_vars = "SCT",
                                               scrna_test_markers = "MAST",
                                               scrna_logfc = 0.585,
                                               scrna_minpct = 0.25,
                                               scrna_label_prefix = "cluster",
                                               pseudocount = 0.083, 
                                               qvalue_cutoff = 0.05,
                                               force_find_markers = FALSE
                                               ){


    # Load Seurat object
    show_msg("Loading {eval(scrna_name)} ..")
    scrna_data <- readRDS(glue::glue("{scrna_folder}/{scrna_name}.rds"))

    # Find or Load Markers
    scrna_markers <- glue::glue("{scrna_folder}/{scrna_name}.{scrna_group_cells_by}.markers")
    Idents(scrna_data) <- scrna_group_cells_by
    if (!file.exists(scrna_markers) | scrna_force_find_markers){

        show_msg("Finding markers for {eval(scrna_name)} ..")
        group_markers <- FindAllMarkers(scrna_data, 
                                          assay = scrna_assay_markers, 
                                          test.use = scrna_test_markers, 
                                          min.pct = scrna_minpct,                      
                                          logfc.threshold = scrna_logfc,
                                          only.pos = F, 
                                          verbose = F)

        write_tsv(group_markers, scrna_markers)

    } else {

        show_msg("Loading markers for {eval(scrna_name)} ..")
        group_markers <- read_tsv(scrna_markers, show_col_types = F)

    }

    # Keep only significant markers
    sig_markers <- group_markers %>% filter(p_val_adj < qvalue_cutoff) %>% pull(gene) %>% unique()
    show_msg("Found {length(sig_markers)} differentially expressed genes at FDR < {eval(qvalue_cutoff)} ..")

    # Also include the most variable genes overall (if they're not included in the above list already)
    var_vars <- VariableFeatures(scrna_data, assay = scrna_assay_vars)
    show_msg("Found {length(var_vars)} variable features ..")

    expressed_genes_of_interest <- union(var_vars,sig_markers)
    show_msg("We will use a total of {length(expressed_genes_of_interest)} expression features ..")


    show_msg("Generating pseudobulk profiles for {eval(scrna_name)} ..")
    get_collapsed_counts <- function(study, label_slot, label_prefix, features=NULL, expression_slot = "RNA"){

        group_name = study[[label_slot]][1,1] %>% as.character()

        if(!is.na(as.numeric(group_name)) && as.numeric(group_name) < 10) group_name = paste0("0",group_name)

        label = paste0(label_prefix,"_",group_name)

        if (expression_slot == "alra"){
            exp_matrix <- study@assays$alra@data %>% as.matrix() %>% t()
        } else {
            exp_matrix <- study@assays$RNA@counts %>% as.matrix() %>% t()
        }

        collapsed <- colSums(exp_matrix) %>% tibble(gene = names(.), counts = .)
        colnames(collapsed) <- c("gene",label)

        if (!is.null(features)) collapsed <- collapsed %>% filter(gene %in% features)

        return(collapsed)

    }


    # Split Seurat objects based on given metadata column
    split_objs <- SplitObject(scrna_data)

    # Sum collapse counts
    bulk_expression <- purrr::map(split_objs, 
                           get_collapsed_counts, 
                           label_slot = scrna_group_cells_by, 
                           label_prefix = scrna_label_prefix)

    bulk_expression <- bulk_expression %>% purrr::reduce(left_join, by="gene")

    # Calculate CPM
    bulk_expression <- bulk_expression %>% column_to_rownames("gene") %>% as.matrix() %>% edgeR::cpm(log = T, prior.count = pseudocount)

    # Zero-out negative expression (if any)
    bulk_expression[bulk_expression < 0] <- 0

    # Remove genes that become all zeros
    bulk_expression <- bulk_expression[which(rowSums(bulk_expression) > 0),]

    # Keep only the relevant features
    expressed_genes_of_interest <- intersect(expressed_genes_of_interest, rownames(bulk_expression))
    bulk_expression <- bulk_expression[expressed_genes_of_interest,]

    # Transpose
    bulk_expression <- bulk_expression %>% t()
    show_msg("Done")
    show_msg("{c('Samples','Features')}: {dim(bulk_expression)}")
    
    return(bulk_expression)

}

