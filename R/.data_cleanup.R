data_cleanup <- function(data) {
  
  # eliminate columns with no measurements
  yy <- colSums(data)
  yy <- which(yy == 0)
  if (length(yy) != 0) {
    data <- data[, -yy]
  }
  
  # keep only variant genes in the data
  # to do that, eliminate the bottom 10% of genes in variance terms
  xx <- sort(apply(data, 1, var), index.return = TRUE)
  z <- ceiling(nrow(data) * .1)
  t <- xx$ix[-seq(1, z)]
  return(t)
}