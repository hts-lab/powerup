#' A test Function
#'
#' This function tests making a function.
#' @param name Your name
#' @keywords test
#' @export
#' @examples
#' mixmap_test("John")
mixmap_test <- function(name = NULL){
    
    if(!is.null(name)){
        
        print(paste0("Hello ",name,"!"))
        
    } else {
        
        print("Hello!")
        
    }
    
}

