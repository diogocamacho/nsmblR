#' Ensemble model (wrapper for nsmblR)
#' 
#' This function will generate an ensemble gene regulatory network model based on different inference algorithms. 
#' 
#' @param data Gene expression data. Matrix is NxM, with genes on the rows and samples on the column.
#' @param gene_names Names for the genes in the data set in the same order as the expression data matrix.
#' @return Returns a ranked data frame in accordance to the \link{\code{consensus}} function.
ensemble_model <- function(data, gene_names) {

  if (missing(data)) stop("Need data.")
  # if (missing(method)) method <- seq(1, 6)
  
  message("---- Network Inference Ensemble Model ----")
  message("")
  
  ui <- .user_inputs()
  
  # message("Data checks...")
  # if (ui[[1]] == 1) {
  #   message("Data is RNA-seq. Computing log2.")
  #   data <- log2(data)
  #   x <- .data_cleanup(data)
  #   data <- data[x, ]
  #   gene_names <- gene_names[x]
  # } else if (ui[[1]] == 2 & ui[[2]] == "N") {
  #   message("Data is microarray and not log2.")
  #   data <- log2(data)
  #   x <- .data_cleanup(data)
  #   data <- data[x, ]
  #   gene_names <- gene_names[x]
  # } else {
  #   x <- .data_cleanup(data)
  #   data <- data[x, ]
  #   gene_names <- gene_names[x]
  # }
  # message("")
  
  message("Processing data matrix...")
  D <- as.matrix(data)
  rownames(D) <- gene_names
  
  message("Calculating mutual information matrix (necessary for CLR, ARACNe, MRNET, and MRNETB)...")
  M <- compute_mi(D)

  message("Inferring networks...")
  m <- seq(1, 7)
  N <- vector(mode = "list", length = length(m))
  for (i in m) {
    if (m[i] == 1 | m[i] == 2 | m[i] == 7) {
      N[[i]] <- infer_network(method = m[i], data = D)
    } else if (m[i] == 3 | m[i] == 4 | m[i] == 5 | m[i] == 6) {
      N[[i]] <- infer_network(method = m[i], data = M)
    }
  }

  message("Merge inference results...")
  mod <- merge_results(network_list = N)
  
  message("Regulatory filtering...")
  reg_mod <- regulatory_filtering(ensemble_df = mod, organism = ui$organism)
  
  message("Filter edges...")
  fed <- edge_filtering(ensemble_df = reg_mod)
  
  message("Edge voting...")
  tt <- edge_voting(ensemble_df = fed)
  
  message("Extracting consensus...")
  cnet <- consensus(vote_tally = tt)
  
  return(cnet)
}