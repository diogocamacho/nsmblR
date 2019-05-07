#' Ensemble model (wrapper for nsmblR)
#' 
#' This function will generate an ensemble gene regulatory network model based on different inference algorithms. 
#' 
#' @param data Gene expression data. Matrix is NxM, with genes on the rows and samples on the column.
#' @param gene_names Names for the genes in the data set in the same order as the expression data matrix.
#' @return Returns a list containing the inferred networks, the filtered networks (based on regulatory interactions) and a ranked data frame in accordance to the \code{\link{consensus}} function.
ensemble_model <- function(data, gene_names) {

  if (missing(data)) stop("Need data.")

  message("---- Network Inference Ensemble Model ----")
  message("")
  
  ui <- user_inputs()
  
  #####
  message("Processing data matrix...")
  D <- as.matrix(data)
  rownames(D) <- NULL
  
  cd <- data_cleanup(data = D, var_thr = 0.1)
  
  if (length(cd$nix_cols) != 0) {
    D <- D[, -cd$nix_cols]
  }
  
  if (length(cd$nix_rows) != 0) {
    D <- D[-cd$nix_rows, ]
    gene_names <- gene_names[-cd$nix_rows]
  }
  
  #####
  message("Calculating mutual information matrix (necessary for CLR, ARACNe, MRNET, and MRNETB)...")
  M <- compute_mi(D)

  #####
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

  #####
  message("Merge inference results...")
  mod <- merge_results(network_list = N)
  
  #####
  message("Regulatory filtering...")
  reg_mod <- regulatory_filtering(ensemble_df = mod, 
                                  organism = ui$organism, 
                                  gene_names = gene_names)

  #####
  message("Edge voting...")
  tt <- edge_voting(ensemble_df = reg_mod)
  
  #####
  message("Extracting consensus...")
  cnet <- consensus(vote_tally = tt)
  
  cnet <- cnet %>%
    dplyr::mutate(., x = replace(x, values = gene_names[as.numeric(x)])) %>%
    dplyr::mutate(., y = replace(y, values = gene_names[as.numeric(y)]))
  
  #####
  message("Writing results...")
  res <- list(inferred_networks = N,
              filtered_networks = reg_mod,
              consensus_network = cnet)
  
  message("DONE.")
  return(res)
}