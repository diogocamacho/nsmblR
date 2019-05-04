user_inputs <- function() {
  
  # organism
  message("Organism:")
  message("[1] E. coli")
  message("[2] Human")
  org <- readline(prompt = "Selection: ")
  if (org == "") {
    stop("Need a valid organism.")
  } else {
    org <- as.integer(org)
    if (org > 2 | org < 1) stop("Invalid organism.")
  }
  
  # data type
  # message("Data type:")
  # message("[1] RNA-seq")
  # message("[2] Microarray")
  # dtype <- readline(prompt = "Selection: ")
  # if (dtype == "") {
  #   message("No data type chosen. Setting to default (log2 normalized microarray data)")
  #   dtype <- 2
  # } else {
  #   dtype <- as.integer(dtype)
  #   if (dtype > 2 | dtype < 1) stop("Invalid data type.")
  # }
  
  
  
  # egde selection method
  # message("Edge selection:")
  # message("[1] Quantile")
  # message("[2] Top edges")
  # message("[3] None")
  # esel <- readline(prompt = "Selection: ")
  # if (esel == "") {
  #   message("No edge selection method chosen. Setting to default (quantile).")
  #   esel <- 1
  # } else {
  #   esel <- as.integer(esel)
  #   if (esel > 3 | esel < 1) stop("Invalid selection.")
  #   if (esel == 1) {
  #     thr <- readline(prompt = "Quantile: ")
  #     thr <- as.numeric(thr)
  #     if (thr < 0 | thr > 1) stop("Invalid quantile.")
  #   } else if (esel == 2) {
  #     thr <- readline(prompt = "Max number of edges: ")
  #     thr <- as.numeric(thr)
  #     if (thr == 0) stop("Cannot select 0 edges.")
  #   } else if (esel == 3) {
  #     thr <- 0
  #   }
  # }
  
  
  # voting method
  # message("Voting method:")
  # message("[1] Majority vote")
  # message("[2] Supermajority")
  # message("[3] Quorum vote")
  # message("[4] Absolute majority")
  # message("[5] Borda count")
  # vmethod <- readline(prompt = "Selection: ")
  # if (vmethod == "") {
  #   message("No voting method selected. Setting to default (Borda count).")
  #   vmethod <- 5
  # } else {
  #   vmethod <- as.integer(vmethod)
  #   if (vmethod < 1 | vmethod > 5) stop("Invalid method.")
  # }
  
  ui <- list(#data_type = dtype, 
             # data_normalization = dnorm, 
             # edge_selection = list(selection_method = esel, threshold = thr),
             organism = org)#,
             #vote_method = vmethod)
  return(ui)
  
}