#' Edge filtering
#' 
#' This function performs a filtering of the inference results based on a quantile calculation of the results. This is currently hard-coded so that the function selects the top 10% edges in each of the inference algorithms.
#' 
#' @param ensemble_df An ensemble set of inference algorithms as defined in the \code{\link{merge_results}} function
#' @return A filtered data frame with edges being quantified (1) or not (0) based on the quantiles for each of the methods.
edge_filtering <- function(ensemble_df) {
  
  fdf <- ensemble_df %>% 
    dplyr::mutate(., cor = replace(cor, which(abs(cor) < quantile(cor, 0.9, na.rm = TRUE)), 0)) %>%
    dplyr::mutate(., clr = replace(clr, which(clr < quantile(clr, 0.9, na.rm = TRUE)), 0)) %>%
    dplyr::mutate(., aracne = replace(aracne, which(aracne < quantile(aracne, 0.9, na.rm = TRUE)), 0)) %>%
    dplyr::mutate(., mrnet = replace(mrnet, which(mrnet < quantile(mrnet, 0.9, na.rm = TRUE)), 0)) %>%
    dplyr::mutate(., mrnetb = replace(mrnetb, which(mrnetb < quantile(mrnetb, 0.9, na.rm = TRUE)), 0)) %>%
    dplyr::mutate(., pcit = replace(pcit, which(pcit < quantile(pcit, 0.9, na.rm = TRUE)), 0)) %>%
    dplyr::mutate(., mutrank = replace(mutrank, which(mutrank < quantile(mutrank, 0.9, na.rm = TRUE)), 0))
  
  fdf[is.na(fdf)] <- 0
  
  return(fdf) # <-- filtered data frame given edge selection method
}