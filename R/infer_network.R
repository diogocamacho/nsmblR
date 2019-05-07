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
#' @param method Integer defining the method to be used. Defaults to using CLR (method 3).
#' @param data Gene expression data frame or mutual information matrix
#' @param quantile_thr Quantile based threshold to threshold inferred networks. Defaults to the top 3% set of edges.
#' @return An inferred network matrix
infer_network <- function(method = c(1, 2, 3, 4, 5, 6, 7), data, quantile_thr) {
  
  if(missing(method)) method <- 3
  if (missing(data)) stop("Need data matrix.")
  if (missing(quantile_thr)) quantile_thr <- 0.97
  
  if (method == 1) {
    message("Method: Spearman correlations")
    net <- corrr::correlate(t(data), method = "spearman", quiet = TRUE) 
    net <- scale(abs(net[, -1]))
    net <- (net + t(net)) / sqrt(2)
    net <- net %>%
      corrr::as_cordf() %>%
      corrr::shave(., upper = TRUE) %>% 
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y) %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else if (method == 2) {
    message("Method: PCIT")
    net <- netbenchmark::pcit.wrap(data = t(data))
    net <- scale(abs(net))
    net <- (net + t(net)) / sqrt(2)
    net <- net %>%
      corrr::as_cordf() %>%
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()  %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y) %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else if (method == 3) {
    message("Method: CLR")
    net <- scale(data)
    net <- (net + t(net)) / sqrt(2)
    net <- net %>% 
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch()  %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y) %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else if (method == 4) {
    message("Method: ARACNe")
    net <- minet::aracne(mim = data, eps = 0.1) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)  %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else if (method == 5) {
    message("Method: MRNET")
    net <- minet::mrnet(data) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)  %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else if (method == 6) {
    message("Method: MRNETB")
    net <- minet::mrnetb(data) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y)))) %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)  %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else if (method == 7) {
    message("Method: MutRank")
    net <- netbenchmark::mutrank.wrap(data = t(data)) %>%
      corrr::as_cordf() %>% 
      corrr::shave(., upper = TRUE) %>%
      corrr::stretch() %>%
      dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
      dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
      dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
      dplyr::arrange(., x, y)  %>% 
      tibble::add_column(., edge = 0) %>%
      dplyr::mutate(., edge = replace(edge, which(abs(r) > quantile(x = abs(r), probs = quantile_thr, na.rm = TRUE)), 1))
  } else {
    stop("Unknown method.")
  }
  return(net)
}