source("~/project/atlasqtl_code/data_simulation/utils.R")

set_pattern <- function(d, p, ind_d0, ind_p0, vec_prob_sh) {
  
  pat <- matrix(FALSE, nrow = p, ncol = d)
  
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
                                     family = "gaussian", pve_per_snp = NULL,
                                     max_tot_pve = NULL, user_seed = NULL) {
  
  if (!is.null(user_seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(user_seed)
  }
  stopifnot(family %in% c("gaussian", "binomial"))
  if (!is.null(pve_per_snp) & !is.null(max_tot_pve))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")
  if (is.null(pve_per_snp) & is.null(max_tot_pve))
    warning(paste("As both pve_per_snp or max_tot_pve were provided as NULL, the ",
                  "pve per SNP was set to its maximum value so that the total ",
                  "pve for the responses are all below 1.", sep = ""))
  
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
  
  if (!is.null(pve_per_snp) & !is.null(max_tot_pve))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")
  
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
                                   pat, vec_maf, pve_per_snp,
                                   max_tot_pve, var_err, chunks_ph = ind_bl, seed = user_seed)
    with(list_eff, {
      phenos <- phenos + snps %*% beta
      if (family == "binomial")
        phenos <- ifelse(phenos > 0, 1, 0)
      list_data <- create_named_list_(phenos, snps, beta,
                                      pat, pve_per_snp)
      class(list_data) <- "sim_data"
      list_data
    })
  })
}

generate_eff_sizes <- function(d, phenos_act, snps_act, pat, vec_maf, pve_per_snp,
                               max_tot_pve, var_err, chunks_ph, seed = seed) {

  set.seed(seed)
  # pve_per_snp average variance explained per snp
  p <- length(vec_maf)
  
  var_snps_act <- apply(snps_act, 2, var)
  var_phenos_act <- apply(phenos_act, 2, var)
  bool_cst <- var_phenos_act == 0
  
  max_per_resp <- max(colSums(pat)) # maximum of number of SNPs assoc. with one protein
  eps <- .Machine$double.eps^0.75
  
  if (is.null(pve_per_snp)) {
    # sets pve_per_snp to the max possible so that the tot_pve for all responses are below 1.
    if (is.null(max_tot_pve)) max_tot_pve <- 1 - eps
    
    pve_per_snp <- max_tot_pve / max_per_resp 
    # divide the total proportion of variance by the number of assoc. proteins (maximum)
  } else {
    
    if (max_per_resp * pve_per_snp > 1 - eps)
      stop(paste("Provided average proportion of variance explained per SNP too ",
                 "high, would lead to a total genetic variance explained above ",
                 "100% for at least one response. \n Setting pve_per_snp < 1 / length(ind_p0) ",
                 "will work for any pattern. \n", sep = ""))
  }
  
  beta <- matrix(0.0, nrow = p, ncol = d)
  ind_d0 <- which(colSums(pat)>0)
  
  beta[, ind_d0] <- sapply(ind_d0, function(k) {
    
    p0_k <- sum(pat[,k]) #number of SNPs assoc. with active protein k
    
    #generate proportion of variation explained by each SNP via Beta distribution
    vec_pve_per_snp <- rbeta(p0_k, shape1 = 2, shape2 = 5) 
    # positively skewed Beta distribution,to give more weight to smaller effect sizes
    
    #scaling it to the upper bound we want
    vec_pve_per_snp <- vec_pve_per_snp / sum(vec_pve_per_snp) * pve_per_snp * p0_k
    
    
    tot_var_expl_k <- pve_per_snp * p0_k * var_err[k] / (1 - pve_per_snp * p0_k)
    
    vec_maf_act <- vec_maf[pat[,k]]
    vec_var_act <- 2 * vec_maf_act * (1 - vec_maf_act)
    
    beta_k <- rep(0.0, p)
    beta_k[pat[,k]] <- sqrt((tot_var_expl_k + var_err[k]) * vec_pve_per_snp / vec_var_act)
    
    # switches signs with probabilty 0.5
    beta_k[pat[,k]] <- sample(c(1, -1), p0_k, replace = TRUE) * beta_k[pat[,k]]
    
    beta_k
  })
  
  create_named_list_(beta, pve_per_snp)
  
}
