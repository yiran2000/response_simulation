library(dplyr)
library(data.table)
require(echoseq)
source("data_simulation_tools/data_generation_utils.R")

#simulate one single hotspot replicating a known hotspot
simulate_SingleHotspot = function(X, Y,
                                  #specify the active SNP
                                  active_SNP = NULL,
                                  hotspot_propensity = 0.1,
                                  # list of active SNPs
                                  active_protein = NULL,
                                  # second shape parameter for the beta distribution controlling the hotspot propensities
                                  sh2 = 20, 
                                  # mean value of beta
                                  mean_beta, 
                                  # user-provided residual correlation matrix
                                  residual_cor_mat = NULL, 
                                  seed){
  if(!is.null(seed)){set.seed(seed = seed)}
  
  X = X_scaled
  Y = Y_scaled
  
  snp_ls = colnames(X)
  protein_ls = colnames(Y)
  
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  
  # whether we should randomly set up active SNP (single) and active protein (multiple)
  if(is.null(active_SNP)){
    active_SNP = sample(snp_ls, 1)
  }
  
  #specify association matrix with BETA
  pat = matrix(nrow = p, ncol = q, 0)
  rownames(pat) = snp_ls
  colnames(pat) = protein_ls
  # pat[active_SNP, active_protein] = 1
  
  if(is.null(active_protein)){
    pat[active_SNP, ] <- sample(c(1,0),length(protein_ls),replace = T,
                                prob = c(hotspot_propensity, 1 - hotspot_propensity))
  }else{
    pat[active_SNP, active_protein] = 1
  }
  

  
  
  # #simulate the BETA values of the active proteins by scaling Beta(2,sh2) to the range of specified mean
  # temp_BETA = rbeta(length(active_protein), shape1 = 2, shape2 = sh2) 
  # BETA = temp_BETA * mean_beta/mean(temp_BETA)
  
  #simulate BETA 
  
  BETA_mat = pat
  row_indices <- match(active_SNP, rownames(pat))
  col_indices <- match(active_protein, colnames(pat))
  BETA_mat[row_indices, col_indices] <- BETA

  
  #replicate the residual_cor_mat
  # R <- Matrix::nearPD(residual_cor_mat, corr = TRUE, do2eigen = TRUE)$mat
  L <- t(chol(residual_cor_mat))
  tZ <- matrix(rnorm(q * n), nrow = q, ncol = n)
  residuals_sim <- t(L %*% tZ) # Gaussian variables
  
  #simulate Y
  # Y_sim = X %*% BETA_mat + residuals_sim
  Y_sim = X %*% BETA_mat + residuals
  
  return(list(
    Y = Y_sim, 
    pat = pat
  ))
  
}


# simulate_withRealGenome(X, p=50, q=30, seed = 1)
