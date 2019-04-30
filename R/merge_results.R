#' Merge results
#' 
#' This function binds all the results from the inference algorithms into a single tidy data frame.
#' 
#' @param network_list A list of inferred networks as defined in the \link{\code{infer_networks}} function.
#' @return A data frame with all the results for the inference models used.
merge_results <- function(network_list) {
  # order interactions
  ordered_list <- lapply(network_list, dplyr::arrange, x, y)
  
  # extract value
  ordered_values <- lapply(ordered_list, dplyr::select, r)
  
  # bind all values
  bound_values <- dplyr::bind_cols(ordered_values)
  colnames(bound_values) <- c("cor", "pcit", "clr", "aracne", "mrnet", "mrnetb", "mutrank")

  # write data frame
  network_dataframe <- dplyr::select(ordered_list[[1]], x, y) %>%
    dplyr::bind_cols(., bound_values)
    
  return(network_dataframe)
}