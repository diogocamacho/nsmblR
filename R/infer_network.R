#' Infer network
#' 
#' This function infers a gene expression network based on pre-defined methods. These are:
#' 1. Spearman correlations
#' 2. Partial correlations with information theory (PCIT)
#' 3. Context Likelihood of Relatedness (CLR)
#' 4. ARACNe
#' 5. MRNET
#' 6. MRNETB
#' 7. Mutual Rank (MutRank)
#' 
#' @param method Integer defining the method to be used
#' @param data Gene expression data frame or mutual information matrix
#' @return An inferred network matrix
infer_network <- function(method = c(1, 2, 3, 4, 5, 6, 7),
                           data) {
  
  if (missing(data)) stop("Need data matrix.")
  
  if (method == 1) {
    # correlations
    message("Method: Spearman correlations")
    net <- corrr::correlate(t(data), method = "spearman",quiet = TRUE) %>%
      corrr::shave(upper = TRUE) %>% 
      corrr::stretch()
  } else if (method == 2) {
    # PCIT
    message("Method: PCIT")
    # y2 <- as.matrix(data)
    # rownames(y2) <- gene_ids
    net <- netbenchmark::pcit.wrap(data = t(data)) %>%
      corrr::as_cordf() %>%
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()  
  } else if (method == 3) {
    # CLR
    message("Method: CLR")
    Z <- scale(data)
    Z <- (Z + t(Z)) / sqrt(2)
    
    net <- corrr::as_cordf(Z) %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()  
  } else if (method == 4) {
    # ARACNE
    message("Method: ARACNe")
    net <- minet::aracne(mim = data, eps = 0.1) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()
  } else if (method == 5) {
    # MRNET
    message("Method: MRNET")
    net <- minet::mrnet(data) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()
  } else if (method == 6) {
    # MRNETB
    message("Method: MRNETB")
    net <- minet::mrnetb(data) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()
  } else if (method == 7) {
    message("Method: MutRank")
    net <- netbenchmark::mutrank.wrap(data = t(data)) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()
  } else {
    stop("Unknown method.")
  }
  
  return(net)
}