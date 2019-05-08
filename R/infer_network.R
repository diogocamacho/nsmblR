#' Infer network
#' 
#' This function infers a gene expression network based on pre-defined methods. 
#' 
#' @details Inference methods:
#' \enumerate{
#'   \item Spearman correlations
#'   \item Partial correlations with information theory (PCIT)
#'   \item Context Likelihood of Relatedness (CLR)
#'   \item ARACNe
#'   \item MRNET
#'   \item MRNETB
#'   \item Mutual Rank (MutRank)
#' }
#' 
#' @param method Integer defining the method to be used. Defaults to using CLR (method 3).
#' @param data Gene expression data frame or mutual information matrix
#' @param quantile_thr Quantile threshold for edge assignment (defaults to 0.95)
#' @return An inferred network matrix
infer_network <- function(method = c(1, 2, 3, 4, 5, 6, 7), data, quantile_thr) {
  
  if (missing(method)) method <- 3
  if (missing(data)) stop("Need data matrix.")
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
  if (method == 1) {
    message("Method: Spearman correlations")
    net <- corrr::correlate(t(data), method = "spearman", quiet = TRUE) %>%
      corrr::shave(., upper = TRUE) %>% 
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., abs(r)) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else if (method == 2) {
    message("Method: PCIT")
    net <- netbenchmark::pcit.wrap(data = t(data)) %>%
      corrr::as_cordf() %>%
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()  %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., abs(r)) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else if (method == 3) {
    message("Method: CLR")
    Z <- scale(data)
    net <- (Z + t(Z)) / sqrt(2)
    net <- corrr::as_cordf(net) %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()  %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., r) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(r >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else if (method == 4) {
    message("Method: ARACNe")
    net <- minet::aracne(mim = data, eps = 0.1) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., r) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(r >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else if (method == 5) {
    message("Method: MRNET")
    net <- minet::mrnet(data) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., r) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(r >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else if (method == 6) {
    message("Method: MRNETB")
    net <- minet::mrnetb(data) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y)))) %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., r) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(r >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else if (method == 7) {
    message("Method: MutRank")
    net <- netbenchmark::mutrank.wrap(data = t(data)) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)
    
    quants <- net %>%
      dplyr::filter(., r != 0) %>% 
      dplyr::select(., r) %>%
      as.matrix %>%
      as.vector %>%
      quantile(probs = quantile_thr, na.rm = TRUE)
    
    net <- net %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(r >= quants), 1)) %>% 
      dplyr::select(., x, y, edge)
    
  } else {
    stop("Unknown method.")
  }
  return(net)
}