#' Data cleanup
#' 
#' Function to clean up the data. The function will eliminate genes with no measurements, genes of low variance, and columns that have no measurement.
#' 
#' @param data Gene expression data matrix
#' @param var_thr Percent variance of genes to be ignored. Eg, setting var_thr = 0.1 will remove the bottom 10% genes in terms of variance
#' @return A list with ids of columns and rows to be removed before continuing analyses.
#' 
#' @examples 
#' Random data:
#' set.seed(123)
#' D <- replicate(expr = ceiling(runif(n = 100, min = 1, max = 1000)), n = 10, simplify = TRUE)
#' cd <- data_cleanup(data = D)
#' 
#' Example data provided:
#' D <- nsmblR::data_matrix
#' cd <- data_cleanup(data = D)
#' 
data_cleanup <- function(data, var_thr) {
  
  if (missing(var_thr)) var_thr <- 0.1
  
  
  # eliminate genes with no measurements
  nix1 <- rowSums(data)
  nix1 <- which(nix1 == 0)

  # eliminate columns with no measurements
  nix2 <- colSums(data)
  nix2 <- which(nix2 == 0)

  # keep only variant genes in the data
  # to do that, eliminate the bottom 10% of genes in variance terms
  nix3 <- sort(apply(data, 1, var), index.return = TRUE)
  z <- ceiling(nrow(data) * var_thr)
  nix3 <- nix3$ix[seq(1, z)]

  # list them
  nix <- list(nix_cols = nix2,
              nix_rows = union(nix1, nix3))
  
  return(nix)
}