clr_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
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
}

spearman_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
  net <- corrr::correlate(t(data), method = "spearman", quiet = TRUE) %>%
    corrr::shave(., upper = TRUE) %>% 
    corrr::stretch() %>%
    dplyr::mutate(., x = replace(x, values = as.numeric(gsub("V", "", x)))) %>% 
    dplyr::mutate(., y = replace(y, values = as.numeric(gsub("V", "", y))))  %>%
    dplyr::mutate(., r = replace(r, list = which(is.na(r)), values = 0)) %>% 
    dplyr::arrange(., x, y)
  
  quants <- net %>%
    dplyr::filter(., r != 0) %>% 
    dplyr::select(., r) %>% 
    abs %>%
    as.matrix %>%
    as.vector %>%
    quantile(probs = quantile_thr, na.rm = TRUE)
  
  net <- net %>% 
    tibble::add_column(., edge = 0) %>%
    dplyr::mutate(., edge = replace(edge, which(abs(r) >= quants), 1)) %>% 
    dplyr::select(., x, y, edge)
}

pcit_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
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
    dplyr::select(., r) %>% 
    abs %>%
    as.matrix %>%
    as.vector %>%
    quantile(probs = quantile_thr, na.rm = TRUE)
  
  net <- net %>% 
    tibble::add_column(., edge = 0) %>%
    dplyr::mutate(., edge = replace(edge, which(abs(r) >= quants), 1)) %>% 
    dplyr::select(., x, y, edge)
}

aracne_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
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
}

mrnet_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
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
}

mrnetb_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
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
}

mutrank_wrapper <- function(data, quantile_thr) {
  
  if (missing(quantile_thr)) quantile_thr <- 0.95
  
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
}