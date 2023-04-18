#' Title
#'
#' @param data A data frame with event names and assigned gates.
#' @param use_gates A vector with gates to be used to create the list.
#'
#' @return A list with gates as names and containing event names.
#' @export
#'
gates_to_list <- function(data, use_gates){

  get_gate_path <- function(gate_name, gate_paths=colnames(data)){

    found <- gate_paths[grep(paste0("\\Q",gate_name,"\\E","$"), gate_paths)]

    if(length(found)>1){
      stop(paste0("Found multiple gates named ", gate_name))
    }
    if(length(found)==0){
      stop(paste0("Found no gate named ", gate_name))
    }

    found

  }

  gate_evid_list <- list()

  for (gate in use_gates) {
    col <- get_gate_path(gate)
    gate_evid_list[[gate]] <- data[data[,col]==TRUE,"EvID_name"]

  }

  gate_evid_list
}
