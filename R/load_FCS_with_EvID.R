#' Load FCS files and add Event ID
#'
#' Adds a new parameter to the FCS file called "EvID" which contains unique event IDs for each event.
#'
#' @param files vector containing FCS file names
#'
#' @return A flowSet
#' @export
#'
#' @importFrom flowCore read.FCS parameters
#' @importFrom BiocGenerics combine
#' @importFrom stats setNames
#' @importFrom methods new as
#'
load_FCS_with_EvID <- function(files){
  flow_list <- list()
  for (fcs in files) {
    print(paste0("Adding Event ID to ", fcs))
    ## Read the original file
    original <- flowCore::read.FCS(fcs, truncate_max_range = FALSE)

    ## Let's create a new parameter as an AnnotatedDataFrame by copying the first parameter from the original flowFrame
    new_p <- flowCore::parameters(original)[1,]

    ## Now, let's change it's name from $P1 to $Px (whatever the next new number is)
    new_p_number <- as.integer(dim(original)[2]+1)
    rownames(new_p) <- c(paste0("$P", new_p_number))

    ## Now, let's combine the original parameter with the new parameter

    allPars <- BiocGenerics::combine(parameters(original), new_p)

    ## Fix the name and description of the newly added parameter, say we want to be calling it cluster_id
    new_p_name <- "EvID"
    allPars@data$name[new_p_number] <- new_p_name
    num.events <- as.integer(dim(original)[1])
    allPars@data$range[new_p_number] <- num.events
    allPars@data$maxRange[new_p_number] <- num.events


    ## Let's get our cluster ID into a single column matrix
    orig_col_names <- dimnames(original@exprs)[[2]]
    cluster_ids <- as.matrix(c(1:num.events), ncol=1)
    new_exprs <- cbind(original@exprs, cluster_ids)
    new_par_col_name <- setNames(new_p_name,
                                 paste0("$P",as.character(new_p_number),"N"))
    dimnames(new_exprs)[[2]] <- c(orig_col_names, new_par_col_name)

    ## Now, let's get all the original keywords and let's add to it
    new_kw <- original@description
    new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P2B"]
    new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
    new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
    new_kw[paste0("$P",as.character(new_p_number),"R")] <- new_kw["$P2R"]
    new_kw[paste0("$P",as.character(new_p_number),"TYPE")] <- "Side_Scatter"
    new_kw[paste0("$P",as.character(new_p_number),"V")] <- new_kw["$P2V"]

    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- 0
    new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- num.events
    new_kw[paste0("P",as.character(new_p_number),"DISPLAY")] <- "LIN"


    ## Now, let's just combine it into a new flowFrame
    new_fcs <- new("flowFrame", exprs=new_exprs, parameters=allPars, description=new_kw)
    out_name <- new_kw["$FIL"][[1]]
    flow_list[[out_name]] <- new_fcs
  }

  return(as(flow_list, "flowSet"))
}
