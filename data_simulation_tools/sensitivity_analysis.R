rm(list = ls())

CORE_DIR <- "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/Projects"

main_dir <- file.path(CORE_DIR, "atlasqtl_addendum/scripts/")
# setwd(main_dir)

RNGkind("L'Ecuyer-CMRG") # to ensure reproducibility when using parallel processes
my_seed <- 123
set.seed(my_seed)

source("data_generation_utils.R")
source("performance_utils.R")

bool_save <- TRUE # whether or not to save the output as RData & pdf files

n_repl <- 50 # number of replicates
n_cpus <- 32 # number of cores

# Generate data #
# ------------- #

require(echoseq)  # package for generation of genotypes and molecular levels; @github/hruffieux/echoseq
require(atlasqtl) # package for joint modelling of QTL hotspots; @github/hruffieux/atlasqtl

n <- 400             # number of samples
q <- 100; q0 <- 50   # number of responses (phenotypes); number of active responses
p <- 1000; p0 <- 20  # number of candidate predictors (SNPs); number of active predictors


max_tot_pve <- 0.05 # upper bound on the proportion of response variance explained by the predictors 
sh2 <- 5           # second shape parameter for the beta distribution controlling the hotspot propensities, 
                   # first shape parameter set to 1, so left-skewed (hotspots mostly small with a few large ones)

real_snps_path <- "" # "../real_snps.rds" # if simulation based on real SNPs, provide a path to and rds file containing a SNP matrix;
                                     # otherwise, e.g., providing the empty string, the SNPs will be simulated

if (file.exists(real_snps_path)) { # use real SNP data (need to be provided by the user)
  
  X <- readRDS(real_snps_path) # must be a matrix of genotypes of size n_real x p_real where n_real >= n and p_real >= p
                               
  stopifnot(is.matrix(X))
  stopifnot(nrow(X) >= n)
  stopifnot(ncol(X) >= p)
  stopifnot(all(X %in% 0:2)) # SNPs must be coded 0-1-2 with no missing values
  
  mess_snps <- "_real_snps"
  
  X <- X[1:n, 1:p, drop = FALSE]
  
} else { # simulate SNPs under Hardy-Weinberg equilibrium and with block-correlation structure 
  
  cor_snps <- c("uncorrelated", "autocorrelated")[2] # whether the SNPs are uncorrelated or autocorrelated by blocks
  
  if (cor_snps == "autocorrelated") {
    
    nb_snps_per_block <- 50
    stopifnot(nb_snps_per_block <= p) # if error, specify a smaller block size
    min_cor <- 0.75; max_cor <- 0.95 # correlation of the Gaussian latent variables used of generating the  
                                     # SNPs drawn uniformly in (min_cor, max_cor) for each block of SNPs
    vec_rho_snps <- runif(ceiling(p / nb_snps_per_block), min = min_cor, max = max_cor) 
    mess_snps <- paste0("_", cor_snps, "_snps_min_", min_cor, "_max_", max_cor)
    
  } else {
    
    cor_snps <- vec_rho_snps <- NULL
    mess_snps <- "_uncorrelated_snps"
    
  }
  
}

cor_phenos <- c("uncorrelated", "equicorrelated")[2] # whether the phenotypes are uncorrelated or equicorrelated by blocks

if (cor_phenos == "equicorrelated") { # residual correlation (blocks of equicorrelated phenotypes)
  
  nb_phenos_per_block <- 10
  stopifnot(nb_phenos_per_block <= q) # if error, specify a smaller block size
  min_cor <- 0; max_cor <- 0.5 # residual correlation of the phenotypes drawn uniformly 
                               # in (min_cor, max_cor) for each block of response
  vec_rho_phenos <- runif(ceiling(q / nb_phenos_per_block), min = min_cor, max = max_cor) 
  mess_phenos <- paste0("_", cor_phenos, "_phenotypes_min_", min_cor, "_max_", max_cor)
  
} else {
  
  cor_phenos <- vec_rho_phenos <- NULL
  mess_phenos <- "_uncorrelated_phenotypes"
  
}

anneal <- c(1, 5, 100) # annealing schedule for atlasqtl 
                       # (geometric schedule, initial temperature, number of temperatures)
bool_save = FALSE
if(bool_save){
  
  res_dir <- paste0(CORE_DIR,"/atlasqtl_addendum/output/sensitivity_hyperparam_n_", n,
                    "_p_", p, "_q_", q, "_p0_", p0, "_q0_", q0, mess_snps, mess_phenos, 
                    "_max_tot_pve_", max_tot_pve, "_sh2_", sh2, 
                    "_anneal_", paste0(anneal, collapse = "_"), "_n_repl_", n_repl, "/")
  dir.create(res_dir)
  
  sink(paste(res_dir, "out.txt", sep=""), append = FALSE, split = TRUE, type = "output")
  sink(file(paste(res_dir, "err.txt", sep=""), open = "wt"), type = "message")
  
}


# SNPs
#
if (file.exists(real_snps_path)) { # prepare the real SNPs in the echoseq format for convenience
  
  rownames(X) <- NULL
  snp_ids <- colnames(X)
  list_X <- NULL
  list_X$snps <- X
  list_X$vec_maf <- apply(X, 2, function(xx) mean(xx) / 2) 
  class(list_X) <- "sim_snps"
  
} else { # simulate the SNPs
  
  list_X <- generate_snps(n = n, p = p, cor_type = cor_snps, vec_rho = vec_rho_snps)
  snp_ids <- colnames(list_X$snps)
  
}

# Association pattern
#
ind_p0 <- sort(sample(1:p, p0, replace = FALSE)) 
ind_q0 <- sort(sample(1:q, q0, replace = FALSE))

vec_prob_sh <- rbeta(p0, shape1 = 1, shape2 = sh2) # hotspot propensities
pat <- set_pattern(q, p, ind_q0, ind_p0, vec_prob_sh)

# Phenotypes
#
list_Y <- generate_phenos(n, q, cor_type = cor_phenos, vec_rho = vec_rho_phenos)

dat <- generate_dependence_pat(list_X, list_Y, pat = pat, max_tot_pve = max_tot_pve)
pat <- dat$pat

list_dat <- list(dat)

if (n_repl > 1) {

  list_dat <- c(list_dat, parallel::mclapply(2:n_repl, function(repl){

    list_Y <- generate_phenos(n, q, cor_type = cor_phenos, vec_rho = vec_rho_phenos)

    generate_dependence_pat(list_X, list_Y, pat = pat, max_tot_pve = max_tot_pve)

  }, mc.cores = n_cpus))

  names(list_dat) <- paste0("repl_", 1:n_repl)
}

require(gplots)
pal <- colorRampPalette(c('yellow','red'))
abs_beta_ind_p0 <- rowSums(sapply(list_dat, function(dats) rowMeans(abs(dats$beta))))[ind_p0] # averaged effect sizes for hotspot
col_pal <- pal(500)[as.numeric(cut(abs_beta_ind_p0, breaks = 500))]

if( bool_save ) {
  pdf(paste0(res_dir, "simulated_hotspot_mean_size.pdf"), width = 7, height = 5, paper='special')
}
display_param(rowSums(pat), 
              ind_p0, 
              col = col_pal, 
              xlab = "SNPs", 
              ylab = "Simulated hotspot mean size", 
              main = "Simulated hotspot pattern")
if (bool_save) {
  dev.off()
}

# Choose hyperparameters #
# ---------------------- #

# Number of predictors associated with each response
# first entry: E_p (prior mean), second entry: V_p (prior variance)
#
list_hyperparam <- list(c(1, 2),
                        c(1, 5), 
                        c(1, 10), 
                        c(2, 5), 
                        c(5, 25), 
                        c(10, 50),
                        c(10, 100)) 

# Run atlasqtl #
# ------------ #

hyper_labels <- sapply(list_hyperparam, function(hyperparam) paste0(hyperparam, collapse = " - "))

list_runs <- parallel::mclapply(1:n_repl, function(repl){
  
  list_atlasqtl <- lapply(list_hyperparam, function(hyperparam) {
    atlasqtl(Y = list_dat[[repl]]$phenos, 
             X = list_dat[[repl]]$snps, 
             p0 = hyperparam, 
             anneal = anneal)$gam_vb }) #
  
  names(list_atlasqtl) <- hyper_labels
  
  list_atlasqtl

}, mc.cores = n_cpus)

# Rearrange results by hyperparameter pairs
#
list_runs <- lapply(hyper_labels, function(hyper_label) {
  list_runs_hyper <- lapply(list_runs, "[[", hyper_label)
  names(list_runs_hyper) <- paste0("repl_", 1:n_repl)
  list_runs_hyper
})
names(list_runs) <- hyper_labels

list_pat <- lapply(1:n_repl, function(repl) pat)

if(bool_save){
  save.image(file = file.path(res_dir, "output.RData"))
  pdf(paste0(res_dir, "ROC.pdf"), width = 6.5, height = 6.5, paper='special')
}
cols <- colorRampPalette(c("black", "grey80"))(length(list_hyperparam))

par(pty="s")
display_all_roc(list_of_lists = list_runs, list_pat = list_pat, 
                title = paste0("ROC curves with 95% CI (p = ", p, ", q = ", q, ")",
                               "\n Selection performance vs hyperparameter choice"),
                title_leg = "E_approx - V_approx", 
                col = cols, leg = hyper_labels)

if(bool_save){
  dev.off()
}


# sessionInfo()
#
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.5
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# Random number generation:
#   RNG:     L'Ecuyer-CMRG 
#  Normal:  Inversion 
#  Sample:  Rejection 
#  
# locale:
# [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ROCR_1.0-7     gplots_3.0.1.1 atlasqtl_0.1.2 echoseq_0.2.3 
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3          mvtnorm_1.0-11      lattice_0.20-38    
#  [4] gtools_3.8.1        PowerTOST_1.4-8     bitops_1.0-6       
#  [7] grid_3.6.1          plyr_1.8.5          KernSmooth_2.23-16 
# [10] TeachingDemos_2.10  cubature_2.0.3      gdata_2.18.0       
# [13] Matrix_1.2-17       tools_3.6.1         RcppEigen_0.3.3.7.0
# [16] parallel_3.6.1      compiler_3.6.1      gsl_2.1-6          
# [19] caTools_1.17.1.2 
