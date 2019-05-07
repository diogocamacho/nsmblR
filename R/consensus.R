#' Consensus network
#' 
#' Based on the tallied votes using the \code{\link{edge_voting}} function, it provides the final consensus network.
#' 
#' @param vote_tally Data frame of tallied votes as defined in the \code{\link{edge_voting}} function
#' @return A  data frame, filtered by majority vote decision (51% or better).
consensus <- function(vote_tally) {
  
  kept_edges <- vote_tally %>%
    dplyr::filter(., majority != 0)

  return(kept_edges)
  
}