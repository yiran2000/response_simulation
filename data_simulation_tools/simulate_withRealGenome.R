library(dplyr)
library(data.table)
require(echoseq)
# source("data_simulation_tools/data_generation_utils.R")
source(file.path(simulation_tool_path, "data_simulation_tools/data_generation_utils.R"))

simulate_withRealGenome = function(X, q = 3000, 
                                   protein_ls = NULL,
                                   active_SNP = NULL,
                                   active_protein = NULL,
                                   active_ratio_p = NULL, 
                                   active_ratio_q = NULL,
                                   residual_cor_mat = NULL,
                                   vec_rho_phenos = NULL,
                                   #rate of missingness in Y (in total, with variation in each column)
                                   missing_ratio = 0, 
                                   m = NULL,
                                   # mean_imputation = F, # whether simulated pheno should be imputed 
                                   # upper bound on the proportion of response variance explained by the predictors 
                                   max_tot_pve = 0.25, 
                                   # second shape parameter for the beta distribution controlling the hotspot propensities
                                   sh2 = 5,
                                   # whether the phenotypes are uncorrelated or equicorrelated by blocks, "uncorrelated" or "equicorrelated"
                                   cor_phenos = "equicorrelated", 
                                   # If phenotypes correlated, size of each correlated block
                                   nb_phenos_per_block =  10, 
                                   # If phenotypes correlated, min and max of correlated strength per block
                                   min_cor = 0, max_cor = 0.5, # residual correlation of the phenotypes drawn uniformly,
                                   seed){
  if(!is.null(seed)){set.seed(seed = seed)}

  n = nrow(X)
  p = ncol(X)
  snp_ls = colnames(X)

  
  if(is.null(active_SNP)){
    p0 = round(active_ratio_p * p)
    ind_p0 = sample(1:p, p0)
  }else{
    p0 = length(active_SNP)
    ind_p0 = which(snp_ls %in% active_SNP)
  }
  
  if(is.null(active_protein)){
    q0 = round(active_ratio_q * q)
    ind_q0 <- sort(sample(1:q, q0, replace = FALSE))
  }else{
    q0 = length(active_protein)
    ind_q0 = sort(which(protein_ls %in% active_protein))
  }
  
  #-----------------------------------------------------------------------------
  # Modify X into our desired format

  list_X <- NULL
  list_X$snps <- as.matrix(X)
  list_X$vec_maf <- apply(X, 2, function(xx) mean(xx) / 2) 
  class(list_X) <- "sim_snps"
  
  rm(X)
  

  #-------------------------------------------------------------------------------
  # Simulate the association pattern 
  # i.e. create a matrix pat storing which predictors are assoc. with which responses
  vec_prob_sh <- rbeta(length(ind_p0), shape1 = 1, shape2 = sh2) # hotspot propensities for each active SNP
  pat <- set_pattern(q, p, ind_q0, ind_p0, vec_prob_sh)
  
  # The way set_pattern works is:
  # 1. first pair up each active SNP with a subset of the active proteins, based on its propensity score
  # 2. Then make sure every active SNP has one protein assoc. and one protein assoc. with one SNP
  #    i.e. if we notice one SNP has no assoc. protein, we assign a random one from the active proteins, same for SNPs

  #-----------------------------------------------------------------------------
  # Simulate Y
  # We first generate Gaussian residuals with specific correlation structures (block-size 10)
  # vec_rho_phenos is the strength of correlation within each block
  if(!is.null(residual_cor_mat)){ # use specified residuals
    

    L <- chol(residual_cor_mat)
    tZ <- matrix(rnorm(q * n), nrow = q, ncol = n)
    residuals_sim <- crossprod(tZ, L) # Gaussian variables
    
    list_Y <- list(
      phenos = residuals_sim,
      var_err = apply(residuals_sim, 2, var)
    )
    class(list_Y) <- "sim_phenos"
    
  }else{
    if (cor_phenos == "equicorrelated" | cor_phenos == "autocorrelated") { # residual correlation (blocks of equicorrelated phenotypes)
      stopifnot(nb_phenos_per_block <= q) # if error, specify a smaller block size
      if(is.null(vec_rho_phenos)){
        vec_rho_phenos <- runif(ceiling(q / nb_phenos_per_block), min = min_cor, max = max_cor) # define the correlation strength rho for each block
      }
      
    } else {
      cor_phenos <- vec_rho_phenos <- NULL
    }
    
    list_Y <- generate_phenos(n, q, cor_type = cor_phenos, vec_rho = vec_rho_phenos, user_seed = seed)
    
  }
  
  # Then we add dependency structure into X and Y:
  # pve for proportion of variance explained by SNPs, here we set a maximum
  dat <- generate_dependence_pat(list_X, list_Y, pat = pat, max_tot_pve = max_tot_pve, m=m, missing_ratio = missing_ratio, user_seed = seed)
  dimnames(dat$beta) <- list(snp_ls, protein_ls)
  dimnames(dat$pat) <- list(snp_ls, protein_ls)
  colnames(dat$phenos) = protein_ls
  
  
  #-----------------------------------------------------------------------------
  # The final X, Y, pat:
  #
  
  # var_X <- apply(X, 2, var)
  # #X <- X[, var_X != 0]
  # pat <- pat[var_X != 0, ]
  
  return(dat)
  
}


# simulate_withRealGenome(X, p=50, q=30, seed = 1)
