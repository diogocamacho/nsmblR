#' Edge voting
#' 
#' Method to define ensemble method. Given a data frame containing inferred relationships between genes, the method will apply a defined voting method to define  a consensus network.
#' Voting methods:
#' 1. Majority vote: edge is considered present if it is inferred in 51% of the cases
#' 2. Super majority vote: edge is considered present if it is inferred in 2/3 of the cases.
#' 3. Quorum vote: this method sits between between the majority vote and the supermajority. It is defined as N/2 + 1, where N is the number of methods used.
#' 4. Absolute majority vote: edge is present in *all* of the methods only.
#' 
#' Note that the function will perform all of these voting procedures on the data.
#' 
#' @param ensemble_df A data frame with the outcomes of all inference methods, as defined in the \code{\link{merge_results}} or in \code{\link{regulatory_filtering}} functions
#' @return A data frame with the outcome of chosen voting method
edge_voting <- function(ensemble_df) {
  
  # voted_edges <- ensemble_df %>% 
  #   dplyr::select(., -x, -y) %>% 
  #   dplyr::mutate(., num_votes = apply(., 1, function(x) length(which(x != 0))))
    
  # majority
  maj <- ensemble_df %>% 
      dplyr::mutate(., majority = vote_count / 7) %>%
      dplyr::mutate(., majority = replace(majority, list = which(majority > 0.51), values = 1)) %>%
      dplyr::mutate(., majority = replace(majority, which(majority != 1), 0))
  
  # super majority
  super_maj <- ensemble_df %>% 
      dplyr::mutate(., super_majority = vote_count / 7) %>%
      dplyr::mutate(., super_majority = replace(super_majority, list = which(super_majority > (2/3)), values = 1)) %>%
      dplyr::mutate(., super_majority = replace(super_majority, which(super_majority != 1), 0))

  # quorum vote
  qvote <- ensemble_df %>% 
    dplyr::mutate(., quorum = vote_count) %>%
      dplyr::mutate(., quorum = replace(quorum, list = which(quorum < ceiling(7 / 2) + 1), values = 0)) %>%
      dplyr::mutate(., quorum = replace(quorum, which(quorum != 0), 1))
    
  # absolute majority  
  abs_maj <- ensemble_df %>% 
      dplyr::mutate(., absolute_majority = vote_count) %>%
      dplyr::mutate(., absolute_majority = replace(absolute_majority, list = which(absolute_majority != 7), values = 0)) %>%
      dplyr::mutate(., absolute_majority = replace(absolute_majority, which(absolute_majority != 0), 1))
  
  # borda voting
  # mod <- ensemble_df %>% 
  #     dplyr::mutate(., rank_cor = rank(1/cor, na.last = TRUE, ties.method = "first")) %>%
  #     dplyr::mutate(., rank_clr = rank(-clr, na.last = TRUE, ties.method = "first")) %>%
  #     dplyr::mutate(., rank_aracne = rank(1/aracne, na.last = TRUE, ties.method = "first")) %>%
  #     dplyr::mutate(., rank_mrnet = rank(1/mrnet, na.last = TRUE, ties.method = "first")) %>%
  #     dplyr::mutate(., rank_mrnetb = rank(1/mrnetb, na.last = TRUE, ties.method = "first")) %>%
  #     dplyr::mutate(., rank_pcit = rank(1/pcit, na.last = TRUE, ties.method = "first")) %>%
  #     dplyr::mutate(., rank_mutrank = rank(1/mutrank, na.last = TRUE, ties.method = "first"))
  # 
  # raw_votes <- cbind(1/mod$cor, 
  #                      1/mod$aracne, 
  #                      -mod$clr, 
  #                      1/mod$mutrank, 
  #                      1/mod$pcit, 
  #                      1/mod$mrnet, 
  #                      1/mod$mrnetb)
  #   
  #   voting_res <- votesys::create_vote(x = t(raw_votes), xtype = 1)
  #   borda_res <- votesys::borda_method(x = voting_res)
  #   borda_ranks <- rank(borda_res$other_info$count_min, ties.method = "first")
  #   
    
    # vote tally
    total_tally <- tibble::tibble(x = ensemble_df$x,
                                  y = ensemble_df$y,
                                  # borda = borda_ranks,
                                  majority = maj$majority,
                                  super_majority = super_maj$super_majority,
                                  absolute_majority = abs_maj$absolute_majority,
                                  quorum = qvote$quorum)
    
    return(total_tally)
    
}    
    
    