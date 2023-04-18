#' Cluster and Run UMAP on Single Cell Experiment
#'
#' @param sce_object A Single Cell Experiment object.
#' @param use_class A string defining which class of markers to use: \code{"type"/"state"}
#' @param seed A numeric specifying which seed to use.
#' @param n A numeric specifying how many cells to subsample; NULL for all cells.
#'
#' @return A Single Cell Experiment object with clustering and UMAP.
#' @importFrom CATALYST runDR cluster
#' @export
#'
cluster_sce <- function(sce_object, use_class, seed=1000, n=5000){
  runDR(cluster(sce_object, features=use_class, seed=seed),
        "UMAP", cells = n, features = use_class, verbose = T)
}
