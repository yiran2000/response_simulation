library(dplyr)
library(data.table)
library(ggplot2)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(patchwork)

source("data_simulation_tools/simulate_withRealGenome.R")

#read in X
geno_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD"
geno = fread(file.path(geno_path, "imputed_chr1_ld1_filtered.raw"))
X = geno[ ,7:ncol(geno)] #randomly subset n number of samples
colnames(X) = sapply(strsplit(colnames(X), "_"), `[`, 1)

################################################################################
# setting up parameters

#-------------------------------------------------------------------------------
#This set of parameters defines the full-scale data similar to our real scenario, but takes a long time to run
param_table = expand.grid(
  #If you just want to see the manhattan plot of the simulated pattern, can change n smaller to, say, 100
  n = 35000, 
  p = 1000,
  q = 3000,
  active_ratio_q = c(0.01, 0.05, 0.2),
  active_ratio_p = 0.02,
  # max_tot_pve = 0.3,
  max_tot_pve = c(0.1, 0.3),
  sh2 = 20
  # X_chunk_size = 300,
  # n_active_chunk = 3
)

#-------------------------------------------------------------------------------
#' Here are the parameters for a smaller dataset to play with.
#' The AURPC and AUROC may be different from the full-scale dataset, 
#' but the relative performance of different vesions of atlasQTL should stay the same
param_table = expand.grid(
  n = 7000,
  p = 200,
  q = 600,
  active_ratio_q = c(0.03, 0.15, 0.6),
  active_ratio_p = 0.1,
  # max_tot_pve = 0.3,
  max_tot_pve = c(0.1, 0.3),
  sh2 = 20
  # X_chunk_size = 60,
  # n_active_chunk = 3
)

################################################################################
# Generate data with the set of parameters defined and plot
dat_ls = lapply( 1:nrow(param_table), function(i) {
  list = simulate_withRealGenome(X, 
                                 n = 100, #I set n = 100 cuz it doesn't matter for this plot and it will be quicker
                                 p=param_table[i, "p"], 
                                 q=param_table[i, "q"], 
                                 missing_ratio = 0,
                                 # X_chunk_size = param_table[i, "X_chunk_size"],
                                 # n_active_chunk = param_table[i, "n_active_chunk"],
                                 active_ratio_q = param_table[i, "active_ratio_q"], 
                                 active_ratio_p = param_table[i, "active_ratio_p"], 
                                 max_tot_pve = param_table[i, "max_tot_pve"],
                                 sh2 = param_table[i, "sh2"],
                                 seed = NULL)
  pat = list$pat
  dat = data.frame(SNP = 1:nrow(pat), 
                   ProteinCount = rowSums(pat),
                   active_ratio_q = param_table[i, "active_ratio_q"],
                   active_ratio_p = param_table[i, "active_ratio_p"],
                   max_tot_pve = param_table[i, "max_tot_pve"],
                   sh2 = param_table[i, "sh2"])
  return(dat) }) %>% rbindlist()

# Make manhattan plot of different scenarios
dat_ls %>% filter(max_tot_pve == 0.1) %>% 
  ggplot(aes(x = SNP, y = ProteinCount))+
  geom_point()+
  facet_wrap(~active_ratio_q+active_ratio_p, labeller = label_both, nrow = 3)+
  theme_bw()


################################################################################
# Generate atlasQTL results for each scenario
res_ls = lapply(1:nrow(param_table), function(i) {
  list = simulate_withRealGenome(X, n = n, p=p, q=q, 
                                 missing_ratio = 0,
                                 # X_chunk_size = param_table[i, "X_chunk_size"],
                                 # n_active_chunk = param_table[i, "n_active_chunk"],
                                 active_ratio_q = param_table[i, "active_ratio_q"], 
                                 active_ratio_p = param_table[i, "active_ratio_p"], 
                                 max_tot_pve = param_table[i, "max_tot_pve"],
                                 sh2 = param_table[i, "sh2"],
                                 seed = NULL)
  X = list$X
  Y = list$Y
  pat = list$pat
  system.time(obj_atlasqtl <- atlasqtl(as.matrix(Y), as.matrix(X),
                                       p0 = c(1, 4),
                                       user_seed = 1, maxit= 50000,
                                       thinned_elbo_eval = F,
                                       batch = "y",
                                       tol = 0.1,
                                       # anneal_tol = c(1, 0.1),
                                       anneal_tol = NULL,
                                       # anneal = c(1, 2, 10),
                                       anneal = NULL,
                                       burn_in = 20,
                                       epsilon_lb = c(2, 1.5, 0.25),
                                       epsilon_it = c(0.1, 50, 0), #k, x_0, m
                                       partial_elbo = F,
                                       eval_perform = T))
  roc <- PRROC::roc.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat2)), curve =TRUE)
  AUROC = roc$auc
  
  pr <- PRROC::pr.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat2)), curve =TRUE)
  AUPRC = pr$auc.integral
  
  return(
    data.frame(active_ratio_q = param_table[i, "active_ratio_q"],
               active_ratio_p = param_table[i, "active_ratio_p"],
             max_tot_pve = param_table[i, "max_tot_pve"],
             AUPRC = AUPRC,
             AUROC = AUROC)
  )
  
  }) %>% rbindlist()



# res_ls %>% group_by(active_ratio_q,active_ratio_p) %>% 
#   mutate(mean_AUPRC = mean(AUPRC), mean_AUROC = mean(AUROC)) %>% 
#   ggplot(aes(x = active_ratio_q, y = AUPRC, color = as.factor(max_tot_pve)))+geom_point()+geom_line()


