library(dplyr)
library(data.table)
require(echoseq)
source("data_simulation_tools/data_generation_utils.R")


simulate_withRealGenome = function(X, n, p, q = 3000, 
                                   # percentage of predictor and response that are actively associated with each other
                                   active_ratio_p = 0.2, active_ratio_q = 0.4,
                                   # X_chunk_size = 100,
                                   # n_active_chunk = 2,
                                   max_q_assoc = 5, # maximum number of proteins to be associated with one protein
                                   #rate of missingness in Y (in total, with variation in each column)
                                   missing_ratio = 0.2, 
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
  
  # X = X[sample(1:nrow(X), n), 1:p] #randomly select n rows and p columns from X
  if(n > nrow(X)|p > nrow(X)){
    X = X[1:n, 1:p]
  }
  p0 = round(active_ratio_p * p)
  q0 = round(active_ratio_q * q)
  
  if (cor_phenos == "equicorrelated") { # residual correlation (blocks of equicorrelated phenotypes)
    
    stopifnot(nb_phenos_per_block <= q) # if error, specify a smaller block size
    vec_rho_phenos <- runif(ceiling(q / nb_phenos_per_block), min = min_cor, max = max_cor) 
    mess_phenos <- paste0("_", cor_phenos, "_phenotypes_min_", min_cor, "_max_", max_cor)
    
  } else {
    
    cor_phenos <- vec_rho_phenos <- NULL
    mess_phenos <- "_uncorrelated_phenotypes"
    
  }
  
  #-----------------------------------------------------------------------------
  # Modify X into our desired format
  snp_ls = colnames(X)
  
  rownames(X) <- NULL
  snp_ids <- colnames(X)
  list_X <- NULL
  list_X$snps <- as.matrix(X)
  list_X$vec_maf <- apply(X, 2, function(xx) mean(xx) / 2) 
  class(list_X) <- "sim_snps"
  
  #-------------------------------------------------------------------------------
  # Simulate the association pattern 
  # i.e. create a matrix pat storing which predictors are assoc. with which responses
  
  # specify the indices of X and Y that are active
  # X
  
  # First divide X into chunk of certain size
  # X_chunk_ls = split(1:p, ceiling(seq_along(1:p) / X_chunk_size))
  # n_X_chunk = length(X_chunk_ls)
  # # Then set random chunks to active according to the porpotion
  # # ind_chunk_active = sample(n_X_chunk, round(n_X_chunk*active_ratio_p), replace = FALSE)
  # ind_chunk_active = sample(n_X_chunk, n_active_chunk, replace = FALSE)
  # # select indexs from chunk
  # ind_p_active_chunks = lapply(ind_chunk_active, function(i) X_chunk_ls[[i]]) %>% unlist

  # Randomly select p0 number of active SNPs
  # ind_p0 = sample(ind_p_active_chunks, p0)
  ind_p0 = sample(1:p, p0)
  
  # ind_p0 <- sort(sample(1:p, p0, replace = FALSE)) 
  ind_q0 <- sort(sample(1:q, q0, replace = FALSE))
  
  vec_prob_sh <- rbeta(length(ind_p0), shape1 = 1, shape2 = sh2) # hotspot propensities for each active SNP
  pat <- set_pattern(q, p, ind_q0, ind_p0, vec_prob_sh)
  # They way set_pattern works is:
  # 1. first pair up each active SNP with a subset of the active proteins, based on its propensity score
  # 2. Then make sure every active SNP has one protein assoc. and one protein assoc. with one SNP
  #    i.e. if we notice one SNP has no assoc. protein, we assign a random one from the active proteins, same for SNPs
  # browser()
  #-----------------------------------------------------------------------------
  # Simulate Y
  # We first generate Gaussian residuals with specific correlation structures (block-size 10)
  # vec_rho_phenos is the strength of correlation within each block
  list_Y <- generate_phenos(n, q, cor_type = cor_phenos, vec_rho = vec_rho_phenos, user_seed = seed)
  
  # Then we add dependency structure into X and Y:
  # pve for proportion of variance explained by SNPs, here we set a maximum
  dat <- generate_dependence_pat(list_X, list_Y, pat = pat, max_tot_pve = max_tot_pve, user_seed = seed)
  pat <- dat$pat
  
  #-----------------------------------------------------------------------------
  # Set missing values
  # 
  Y = dat$phenos
  if(missing_ratio > 0){
    num_entries <- round(missing_ratio * length(Y))
    indices <- sample(seq_along(Y), num_entries)
    Y[indices] <- NA
  }
  
  #-----------------------------------------------------------------------------
  # The final X, Y, pat:
  #
  X = as.matrix(dat$snps)
  Y = as.matrix(Y)
  
  var_X <- apply(X, 2, var)
  #X <- X[, var_X != 0]
  pat <- pat[var_X != 0, ]
  
  return(list(
    X = X, 
    Y = Y, 
    pat = pat
  ))
  
}


# simulate_withRealGenome(X, p=50, q=30, seed = 1)
