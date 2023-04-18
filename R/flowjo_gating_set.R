#' Apply FlowJo Gating on FlowSet
#'
#' Uses a FlowJo workspace file to apply gating, transformations, etc. on a flowSet.
#'
#' @param flowset A flowSet containing flowFrames.
#' @param flowjo_workspace FlowJo workspace file.
#'
#' @return A GatingSet containing GatingHierarchies.
#' @importFrom flowWorkspace flowSet_to_cytoset gh_apply_to_cs
#' @importFrom CytoML open_flowjo_xml flowjo_to_gatingset
#' @export
#'
get_gating_set <- function(flowset,flowjo_workspace){
  cs <- flowSet_to_cytoset(flowset)
  ws <- open_flowjo_xml(flowjo_workspace)
  gs <- flowjo_to_gatingset(ws, name=1, execute = F) #parsing workspace without loading fcs data

  gh_apply_to_cs(gs[[1]], cs, compensation_source = "template") #applying gating, transformations, etc. to cytoset and making a gatingset
}
