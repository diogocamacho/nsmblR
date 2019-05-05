#' Merge results
#' 
#' This function binds all the results from the inference algorithms into a single tidy data frame.
#' 
#' @param network_list A list of inferred networks as defined in the \code{\link{infer_network}} function.
#' @return A data frame with all the results for the inference models used.
merge_results <- function(network_list) {
  
  # edge calls
  ordered_values <- lapply(network_list, dplyr::select, edge)
  
  # bind all values
  bound_values <- dplyr::bind_cols(ordered_values) %>% dplyr::mutate(tally = rowSums(.))
  colnames(bound_values) <- c("cor", "pcit", "clr", "aracne", "mrnet", "mrnetb", "mutrank", "vote_count")

  # write data frame
  network_dataframe <- dplyr::select(network_list[[1]], x, y) %>%
    dplyr::bind_cols(., bound_values) %>% 
    dplyr::filter(., vote_count != 0)
    
  return(network_dataframe)
}