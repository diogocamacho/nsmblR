#' Regulators filtering
#' 
#' This function will filter the inference results to only those between transcriptional regulators and protein coding genes. It assumes that gene-gene interactions happen due to the effect of a transcriptinal regulator. This function should be run independent of the outcome of edge filtering done with the \code{\link{edge_filtering}} function.
#'
#' @param ensemble_df A data frame from \code{\link{merge_results}}
#' @param organism Organism being analyzed (E. coli or Human)
#' @param gene_names Set of gene names in the data set (ordered as in the expression data matrix)
#' @return A filtered data frame 
regulatory_filtering <- function(ensemble_df, organism, gene_names) {
  
  if(organism == 1) {
    load("data/regulondb.RData")
    regs <- union(regdb$tf_predictions$gene_symbol, regdb$srnas$gene_symbol)
  } else {
    load("data/ttrust.RData")
    regs <- ttrust_data$tf
  }
  
  b1 <- gene_names[ensemble_df$x]
  b2 <- gene_names[ensemble_df$y]
  
  a1 <- which(b1 %in% regs)
  a2 <- which(b2 %in% regs)
  
  tf_df <- ensemble_df[union(a1, a2), ]
  
  return(tf_df)
} 