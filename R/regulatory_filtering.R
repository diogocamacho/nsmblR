#' Regulators filtering
#' 
#' This function will filter the inference results to only those between transcriptional regulators and protein coding genes. It assumes that gene-gene interactions happen due to the effect of a transcriptinal regulator. This function should be run independent of the outcome of edge filtering done with the \link{\code{edge_filtering}} function.
#'
#' @param ensemble_df
#' @return A filtered data frame 
regulatory_filtering <- function(ensemble_df, organism) {
  
  if(organism == 1) {
    load("data/regulondb.RData")
    regs <- union(regdb$tf_predictions$gene_symbol, regdb$srnas$gene_symbol)
  } else {
    load("data/ttrust.RData")
    regs <- ttrust_data$tf
  }
  
  a1 <- which(ensemble_df$x %in% regs)
  a2 <- which(ensemble_df$y %in% regs)
  
  tf_df <- ensemble_df[union(a1, a2), ]
  
  return(tf_df)
} 