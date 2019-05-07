#' Edge filtering
#' 
#' This function performs a filtering of the inference results based on a quantile calculation of the results. This is currently hard-coded so that the function selects the top 10% edges in each of the inference algorithms.
#' 
#' @param inferred_network An inferred networks as defined in the \code{\link{infer_network}} function
#' @param quantile_thr Quantile based threshold to threshold inferred networks. Defaults to the top 10% set of edges.
#' @return A filtered data frame with edges being quantified (1) or not (0) based on the quantiles for each of the methods.
edge_filtering <- function(inferred_network, quantile_thr) {
  
  quants <- inferred_network %>%
    dplyr::filter(., r != 0) %>% 
    dplyr::select(., r) %>%
    as.matrix %>%
    as.vector %>%
    quantile(probs = quantile_thr, na.rm = TRUE)
  
  net <- inferred_network %>% 
    tibble::add_column(., edge = 0) %>%
    dplyr::mutate(., edge = replace(edge, which(r >= quants), 1))
  return(net)
}