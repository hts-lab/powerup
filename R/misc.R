#' Display message
#'
#' This function prints a message and calls flush.console() to refresh output
#' @param x Message to be shown, will be used inside glue() call.
#' @keywords message
#' @export
#' @examples
#' show_msg("The number of cells is {cell_count}")
show_msg <- function(x){
    
    cat(glue::glue(x),sep = "\n")
    flush.console()
    
}

