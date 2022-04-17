#' Display message
#'
#' This function prints a message and calls flush.console() to refresh output
#' @param x Message to be shown, will be used inside glue() call.
#' @keywords message
#' @import glue
#' @export
#' @examples
#' show_msg("The number of cells is {cell_count}")
show_msg <- function(x){
    
    cat(glue::glue(x, .envir = parent.frame()),sep = "\n")
    flush.console()
    
}



#' Some random new function
#'
#' This function prints "hello" to test something.
#' @param x Name of person to say hello to.
#' @keywords hello
#' @export
#' @examples
#' say_hello("Mush")
say_hello <- function(x){
  
  cat(paste("Hello ",x),sep = "\n")
  flush.console()
  
}




#' Extract ggplot legend
#'
#' This function extracts the ggplot legend to plot separately.
#' @param ggplot_object A ggplot object.
#' @keywords ggplot legend
#' @export
#' @examples
#' get_legend(my_ggplot)
get_legend <- function(ggplot_object){
    
   tmp <- ggplot_gtable(ggplot_build(ggplot_object))
   legend_id <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[legend_id]]
   return(legend)
                      
}
  

#' Get a palette of specified size
#'
#' This function creates a palette of specified size using equal distant hues.
#' @param n The size of the palette.
#' @keywords hues
#' @export
#' @examples
#' gg_color_hue(5)                           
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
  


#' Convert arxspan ID to cell line name
#'
#' This function converts an arxpan ID to a cell line name.
#' @param arxspan The arxspan ID (aka DepMap_ID).
#' @param sample_info The dataset containing the cell line name. Rownames must be arxspan IDs, and cell line name must be in column 'CCLE_Name'.
#' @keywords arxspan
#' @export
#' @examples
#' get_cell_line_name("ACH-000001",sample.info)    
get_cell_line_name <- function(arxspan, sample_info){
    
    return(sample_info[arxspan,"CCLE_Name"])
    
}



#' Changes the value column name of a dataframe
#'
#' This function converts an arxpan ID to a cell line name.
#' @param df The dataframe.
#' @param new_colname The new name.
#' @keywords colname
#' @export
change_colname <- function(df, new_colname){
  df %>% rename(!!new_colname := value)
}


#' Double-digit converter
#'
#' This function converts a single digit number to double digits.
#' @export
get_double_digit_number <- function(x){
  
  suppressWarnings(
    if(!is.na(as.numeric(as.character(x))) && as.numeric(as.character(x)) < 10) return(paste0("0",x))
  )
  
  return(as.character(x))
  
}

#' Vectorized double-digit converter
#'
#' This function converts a single digit number to double digits.
#' @export
Get_Double_Digit_Number <- Vectorize(get_double_digit_number)




#' Get demo samples
#'
#' This is a helper function to automatically pick cell lines with opposite predictions
#' @export
get_demo_samples <- function(model, samples = NULL, lineage = NULL, model_data = NULL){
  
  selection <- model$predictions
  acceptable_names <- names(model$predictions)
  
  if(!is.null(lineage) && !is.null(model_data)){
    data <- get_original_data(model, model_data) %>% filter(get(paste0("Lin_",lineage)) == 1)
    acceptable_names <- rownames(data)
  }
  
  if(!is.null(samples) && length(samples) > 0){
    return(intersect(samples, acceptable_names))
  }
  
  selection <- selection[acceptable_names]
  if(is.null(samples) || length(samples) == 0){
    sample_top <- selection[which(selection == max(selection))] %>% names()
    
    sample_bottom <- selection[which(selection == min(selection))] %>% names() 
    return(c(sample_top,sample_bottom))
  } else {
    return(intersect(samples, acceptable_names))
  }
  
}


