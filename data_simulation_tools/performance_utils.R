source("utils.R")

compute_roc_list <- function(list_ppi, list_pat) {
  
  list_vec_rank <- lapply(list_ppi, function(ppi) as.numeric(as.factor(ppi)))
  list_vec_pat <- lapply(list_pat, function(pat) as.numeric(pat))
  
  require(ROCR)
  pred <- prediction(list_vec_rank , list_vec_pat)
  performance(pred, measure = "tpr", x.measure = "fpr")
  
}

display_all_roc <- function(list_of_lists, list_pat, title = "ROC curves", 
                            title_leg = "", cols = NULL, leg = NULL) {
  
  ll <- length(list_of_lists)
  if (is.null(cols))
    cols <- rainbow(ll)
  
  par(pty="s")
  
  perf <- compute_roc_list(list_of_lists[[1]], list_pat)
  ROCR::plot(perf, main = title,
             col = cols[1], type = "l", lwd = 2, 
             avg = "vertical", spread.estimate = "stderror", spread.scale = 2)
  abline(0,1)
  
  if (ll > 1) {
    for (i in 2:ll) {
      
      perf <- compute_roc_list(list_of_lists[[i]], list_pat)
      ROCR::plot(perf, 
                 col = cols[i], type = "l", lty = 1, lwd = 2, add = T, 
                 avg="vertical", spread.estimate = "stderror", spread.scale = 2)
      
    }
    
  }
  
  if (!is.null(leg)) {
    legend("bottomright", leg, title = title_leg, col = cols, cex = 0.9, 
           lty = 1, lwd = 2, bty = "n")
  }
  
} 

display_param <- function(vec, ind, xlab = "", ylab = "", main = "", cols = NULL) {
  
  if(is.null(cols))
    cols <- "red"
  
  par(mfrow = c(1, 1))
  plot(vec, xlab = xlab, ylab = ylab, main = main)
  points(ind, vec[ind], col = cols, pch = 19)
  points(ind, vec[ind], pch = 21)
  
}