#' Consensus network
#' 
#' Based on the tallied votes using the \link{\code{edge_voting}} function, it provides the final consensus network.
#' 
#' @param vote_tally Data frame of tallied votes as defined in the \link{\code{edge_voting}} function
#' @return A  data frame, filtered by super majority vote decision.
consensus <- function(vote_tally) {
  
  kept_edges <- vote_tally %>%
    dplyr::filter(., super_majority != 0)

  return(kept_edges)
  
}