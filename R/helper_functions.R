#' List all .FCS files in a given path
#'
#' @param fcs_path A string with the path to the .FCS files.
#'
#' @return A vector with the .FCS file names.
#' @export
#'
list_fcs <- function(fcs_path){
  list.files(path=fcs_path, pattern=".fcs$",full.names = TRUE)
}


#' Convert to Factor and Order Levels According to Template
#'
#' @param vect A vector or factor.
#' @param level_order A vector with given order of the factor levels.
#'
#' @return A factor with the specified order of levels.
#' @export
#'
relevel_factor <- function(vect, level_order){
  return(factor(as.factor(vect), levels=level_order))
}
