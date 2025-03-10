# source("data_simulation_tools/utils.R")
source(file.path(simulation_tool_path, "data_simulation_tools/utils.R"))

set_pattern <- function(d, p, ind_d0, ind_p0, vec_prob_sh) {
  
  pat <- matrix(FALSE, nrow = p, ncol = d)
  # browser()
  for(j in 1:length(ind_p0)) {
    
    ind_j <- ind_p0[j]
    pat[ind_j, ind_d0] <- sample(c(TRUE, FALSE),
                                 length(ind_d0), replace = TRUE,
                                 prob = c(vec_prob_sh[j], 1 - vec_prob_sh[j]))
    
    if (all(!pat[ind_j, ind_d0]))
      pat[ind_j, ind_d0][sample(1:length(ind_d0), 1)] <- TRUE 
    # if active SNP j has no associations, then select a random protein and pair it up
  }
  
  for(ind_k in ind_d0) {
    if (all(!pat[ind_p0, ind_k]))
      # each active response must be associated with at least one predictor
      pat[ind_p0, ind_k][sample(1:length(ind_p0), 1)] <- TRUE
  }
  pat
}

generate_dependence_pat <- function (list_snps, list_phenos, pat,
                                     family = "gaussian", m = NULL,
                                     missing_ratio = 0,
                                     max_tot_pve = NULL, user_seed = NULL) {
  
  if (!is.null(user_seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(user_seed)
  }
  stopifnot(family %in% c("gaussian", "binomial"))
  if ((!is.null(m) & !is.null(max_tot_pve)) | (is.null(m) & is.null(max_tot_pve)))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")
  if (!inherits(list_snps, "sim_snps"))
    stop(paste("The provided list_snps must be an object of class ``sim_snps''. \n",
               "*** You must either use the function generate_snps to simulate snps ",
               "under Hardy-Weinberg equilibrium or the function replicate_real_snps ",
               "to simulate SNPs from real SNP data, by replicating their minor ",
               "allele frequencies and linkage disequilibrium structure. ***",
               sep = ""))
  
  if (!inherits(list_phenos, "sim_phenos"))
    stop(paste("The provided list_phenos must be an object of class ``sim_phenos''. \n",
               "*** You must either use the function generate_phenos to simulate ",
               "phenotypes from (possibly correlated) Gaussian variables or the ",
               "function replicate_real_phenos to simulate phenotypes from real ",
               "phenotypic data, by replicating their correlation structure. ***",
               sep = ""))

  with(c(list_snps, list_phenos), {
    d <- ncol(phenos)
    n <- nrow(snps)
    p <- ncol(snps)
    
    
    if (n != nrow(phenos))
      stop("The numbers of observations used for list_snps and for list_phenos do not match.")
    if (n == 1)
      stop("The number of observations must be greater than 1 in order to generate associations.")
    phenos_act <- phenos[, colSums(pat)>0, drop = FALSE]
    snps_act <- snps[, rowSums(pat)>0, drop = FALSE]

    list_eff <- generate_eff_sizes(d, phenos_act, snps_act,
                                   pat, vec_maf, 
                                   max_tot_pve=max_tot_pve, m=m, 
                                   var_err, seed = user_seed)
    

    with(list_eff, {
      phenos <- phenos + snps %*% beta
      if (family == "binomial")
        phenos <- ifelse(phenos > 0, 1, 0)
      # Set missing values
      if(missing_ratio > 0){
        num_entries <- round(missing_ratio * length(phenos))
        indices <- sample(seq_along(Y), num_entries)
        Y[indices] <- NA
      }
      
      list_data <- create_named_list_(phenos, snps, beta,
                                      pat)
      class(list_data) <- "sim_data"
      list_data
    })
  })
}

generate_eff_sizes <- function(d, phenos_act, snps_act, pat, vec_maf,
                               max_tot_pve = NULL, m = NULL, var_err, seed = seed) {
  
  set.seed(seed)
  # pve_per_snp average variance explained per snp
  p <- length(vec_maf)
  
  var_snps_act <- apply(snps_act, 2, var)
  var_phenos_act <- apply(phenos_act, 2, var)
  bool_cst <- var_phenos_act == 0

  beta <- matrix(0.0, nrow = p, ncol = d)
  ind_d0 <- which(colSums(pat)>0)

  if(!is.null(max_tot_pve)){
    max_per_resp <- max(colSums(pat)) # maximum of number of SNPs assoc. with one protein
    pve_per_snp <- max_tot_pve / max_per_resp 
  }else{
    vec_pve_k = rep(0.0, ncol(pat))
    vec_pve_k[ind_d0] = rbeta(length(ind_d0), 1, (1-m)/m) #m is the desired mean for this BETA distibution
  }
  
  # suppose we have only one active SNP
  if(sum(rowSums(pat)>0) == 1){

    active_snp_ind = which(rowSums(pat)>0)
  
    maf = vec_maf[active_snp_ind]
    var_x <- 2 * maf * (1 - maf)
    
    var_y = var_err/(1-vec_pve_k)
    beta[active_snp_ind, ] = sqrt(var_y * vec_pve_k / var_x)

    beta[active_snp_ind, ind_d0] <- sample(c(1, -1), length(ind_d0), replace = TRUE) * beta[active_snp_ind, ind_d0]

  }else{
    # browser()
    beta[, ind_d0] <- sapply(ind_d0, function(k) {
      
      pve_k = vec_pve_k[k]
      
      p0_k <- sum(pat[,k]) #number of SNPs assoc. with active protein k
      
      #generate proportion of variation explained by each SNP via Beta distribution
      vec_pve_per_snp <- rbeta(p0_k, shape1 = 2, shape2 = 5) 
      # positively skewed Beta distribution,to give more weight to smaller effect sizes
      
      
      #scaling vec_pve_per_snp
      if(!is.null(max_tot_pve)){    
        #case 1: scale it to have mean = max_tot_pve/maxmum number of active SNPs
        vec_pve_per_snp <- vec_pve_per_snp / sum(vec_pve_per_snp) * pve_per_snp * p0_k
        # tot_var_expl_k <- pve_per_snp * p0_k * var_err[k] / (1 - pve_per_snp * p0_k)
        tot_var_k <- var_err[k]/(1- pve_per_snp * p0_k)
      }else{     
        #case 2: scale it to the varied ht^2
        vec_pve_per_snp <- vec_pve_per_snp / sum(vec_pve_per_snp) * pve_k
        tot_var_k <- var_err[k]/(1-pve_k)
      }
      
      
      vec_maf_act <- vec_maf[pat[,k]]
      vec_var_act <- 2 * vec_maf_act * (1 - vec_maf_act)
      
      beta_k <- rep(0.0, p)
      # beta_k[pat[,k]] <- sqrt((tot_var_expl_k + var_err[k]) * vec_pve_per_snp / vec_var_act)
      beta_k[pat[,k]] <- sqrt( tot_var_k * vec_pve_per_snp / vec_var_act)
      
      
      # switches signs with probabilty 0.5
      beta_k[pat[,k]] <- sample(c(1, -1), p0_k, replace = TRUE) * beta_k[pat[,k]]
      
      beta_k
    })
  }
  

  create_named_list_(beta)
  
}

