#' Designate Non-Overlapping Gates to Column Data of a Single Cell Experiment
#'
#' @param sce_obj A Single Cell Experiment Object.
#' @param gate_list A list with gates as names and containing event names.
#' @param gate_levels A vector containing the order of the gates.
#' @param col_name A string defining the column name in the colData of the Single Cell Experiment object.
#'
#' @return A Single Cell Experiment object with the added column containing the designated gates.
#' @importFrom SingleCellExperiment colData
#' @export
#'
unique_gate_designation <- function(sce_obj, gate_list, gate_levels, col_name){

  desig_df <- data.frame(evid=rownames(colData(sce_obj)), desig="other")
  colnames(desig_df)[2] <- col_name

  if (! all(gate_levels %in% names(gate_list))){
    stop("Gate_levels are not in gate_list.")
  }

  for (gate_name in rev(gate_levels)) {
    current_evid <- gate_list[[gate_name]]
    desig_df[desig_df$evid %in% current_evid,col_name] <- gate_name
  }

  desig_df[,col_name] <- relevel_factor(desig_df[,col_name], c(gate_levels,"other"))

  colData(sce_obj)[,col_name] <- desig_df[,col_name]

  sce_obj
}
