#' Assign Downstream Gates to Events from Single Cell Experiment Object
#'
#' @param sce_object A Single Cell Experiment object.
#' @param gatingset A GatingSet with GatingHierarchies.
#' @param use_gate A string defining which gate to use.
#'
#' @return A data frame with event names and assigned gates.
#' @importFrom flowCore exprs
#' @importFrom flowWorkspace gh_pop_get_data gh_pop_get_descendants sampleNames
#' @importFrom purrr map map2
#' @importFrom stats setNames
#' @importFrom SingleCellExperiment colData
#' @export
#'
assign_downstream_gates <- function(sce_object, gatingset, use_gate){
  create_gated_event_list <- function(gating_hierarchy, pregate=use_gate){

    get_EvID <- function(gate_name){
      exprs(gh_pop_get_data(gating_hierarchy, gate_name))[,"EvID"]
    }

    gate_paths <- unlist(map2(gating_hierarchy, pregate, gh_pop_get_descendants))

    gated_evid <- map(gate_paths, get_EvID)

    gated_evid <- map2(paste0(sampleNames(gating_hierarchy),"_"), gated_evid, paste0)
    names(gated_evid) <- gate_paths

    gated_evid
  }

  gated_event_list_all <- map(gatingset, create_gated_event_list) #Takes time

  keys <- unique(unlist(map(gated_event_list_all, names)))
  concat_lists <- function(list_of_lists){
    setNames(do.call(mapply, c(FUN=c, lapply(list_of_lists, `[`, keys))), keys)
  }

  if (length(keys) == 1){
    gated_event_list <- list(unlist(gated_event_list_all, use.names = FALSE))
    names(gated_event_list) <- keys
  } else{
    gated_event_list <- concat_lists(gated_event_list_all)
  }

  gated_evid <- rownames(colData(sce_object))
  gated_event_df <- as.data.frame(matrix(nrow = length(gated_evid),
                                         ncol=length(gated_event_list)))
  colnames(gated_event_df) <- names(gated_event_list)

  for (gate in colnames(gated_event_df)) {
    gated_event_df[,gate] <- gated_evid %in% gated_event_list[[gate]]
  }
  gated_event_df$EvID_name <- gated_evid

  gated_event_df
}
