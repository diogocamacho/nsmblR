#' Compute mutual information matrix
#' 
#' This function will use the `buid.mim` function in the `minet` package to compute the mutal information matrix in the gene space.
#' @param data  Data matrix. This is an NxM matrix, where genes are on the rows and samples on the columns
#' @return Returns the mutual information matrix for the data.
#' 
#' @examples
#' M <- nsmblR::data_matrix
#' compute_mi(M)
compute_mi <- function(data) {
  M <- minet::build.mim(dataset = t(data), estimator = "spearman")
  return(M)
}