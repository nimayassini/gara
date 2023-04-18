#' Apply Gate on GatingHierarchy
#'
#' @param gating_hierarchy A GatingHierarchy.
#' @param gate_name A string with the name of the gate to be used.
#'
#' @return A data frame with the gated data.
#' @importFrom flowWorkspace gh_pop_get_data sampleNames
#' @importFrom flowCore exprs
#'
.gate_data <- function(gating_hierarchy, gate_name){

  out_table <- as.data.frame(exprs(gh_pop_get_data(gating_hierarchy,gate_name)))
  out_table$cell_gate <- gate_name
  out_table$sample_id <- sampleNames(gating_hierarchy)

  out_table
}






#' Asinh Transformation of Expression Data
#'
#' @param xpr_data Expression data.
#' @param channels Paramters to have the transformation applied to.
#' @param cofactor Cofactor to use for transformation.
#'
#' @return Transformed expression data.
#'
.transform_xpr <- function(xpr_data,
                          channels,
                          cofactor=150){

  if (length(cofactor) != length(channels) && length(cofactor) != 1){
    print("Cofactors given do not match length of channels.")
    print("Using cofactor 150 for all channels given!")
    cofactor = 150
  }

  xpr_data[,channels] <- asinh(xpr_data[,channels]/cofactor)

  xpr_data
}







#' Create Single Cell Experiment Object from FlowSet
#'
#' @param flowset A flowSet object containing flowFrames.
#' @param use_gate A vector containing names of gates to use.
#' @param gatingset A GatingSet object to be used to derive Events in defined gates.
#' @param sample_vector A vector with sample names.
#' @param sample_levels A vector defining the levels of the sample conditions.
#' @param panel_df A data frame with the panel information.
#'
#' @return A Single Cell Experiment object with the gated FCS files.
#' @importFrom purrr map map2
#' @importFrom flowCore fsApply keyword compensate
#' @importFrom flowWorkspace gh_get_compensations
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter
#' @export
#'
#'
sce_from_flowset <- function(flowset, use_gate, gatingset,
                             sample_vector, sample_levels, panel_df){

  pass_gs_to_gate <- function(gating_set){
    base::do.call("rbind", map2(gating_set, use_gate, .gate_data))
  }

  message("Generating Gated DataFrame.")
  gated_df <- as.data.frame(base::do.call("rbind", map(gatingset, pass_gs_to_gate)))

  sample_names <- unlist(fsApply(flowset, keyword, "$FIL"))
  names(sample_names) <- c()

  message("Generating list of Event IDs.")
  evid_list <- map2(list(gated_df), sample_names,
                    function(gated_data_frame, sample_name){
                      filter(gated_data_frame, sample_id==sample_name)$EvID
                    })
  message("Preparing Gated flowSet.")
  comp <- gh_get_compensations(gatingset[[1]])@spillover
  rownames(comp) <- sub("Comp-","",rownames(comp))
  colnames(comp) <- sub("Comp-","",colnames(comp))
  fs_gated <- map2(compensate(flowset, comp), evid_list,
                   function(flowframe, evid) {flowframe[evid,]})

  xpr <- do.call("rbind", map(fs_gated, exprs))
  xpr_t <- .transform_xpr(xpr, channels = panel_df[panel_df$class!="none","channel"])
  colnames(xpr_t) <- NULL
  xpr_t <- t(as.matrix(xpr_t))

  # construct colData
  message("Constructing colData.")
  sample_id <- as.factor(gated_df$sample_id)
  cell_gate <- as.factor(gated_df$cell_gate)
  event_id <- as.factor(gated_df$EvID)

  condition <- character()

  for (i in 1:length(sample_names)) {
    cond_index <- which(sample_id==sample_names[i])
    condition[cond_index] <- sample_vector[i]
  }
  condition <- relevel_factor(condition, sample_levels)

  cd <- DataFrame(condition,sample_id,cell_gate,event_id,
                  row.names = paste0(sample_id,"_",event_id))


  print(as_tibble(unique(cd[,c("condition","sample_id")])))


  # construct rowData
  message("Constructing rowData.")
  panel_rd <- panel_df
  panel_rd[is.na(panel_df$antigen), "antigen"] <- panel_rd[is.na(panel_df$antigen), "channel"]
  rd <- DataFrame(
    row.names = panel_rd$antigen, channel_name = panel_rd$channel,
    marker_name = panel_rd$antigen, marker_class = panel_rd$class)


  gate_meta <- data.frame(condition=as.factor(sample_vector),
                          sample_id=as.factor(sample_names))

  #make Single Cell Experiment
  SingleCellExperiment(
    assays = list(exprs = xpr_t),
    rowData = rd, colData = cd,
    metadata = list(experiment_info = gate_meta))


}
