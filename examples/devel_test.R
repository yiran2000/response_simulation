library(profvis)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(tictoc)
library(PRROC)
library(echoseq)
library(peakRAM)


simulation_tool_path = "~/project/response_simulation"
list = readRDS(file.path("/rds/user/yl2021/hpc-work/hotspot_sim", "sim_list_3.rds"))
X_real = as.matrix(list$X)


source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source(file.path(simulation_tool_path, "demo_tools/res_demo_helpers.R"))


rm(list)
gc()

p = 200
q = 1000

list_sim = simulate_withRealGenome(
  X = X_real[,1:p], 
  q = q,
  protein_ls = paste0("protein_", 1:q),
  active_protein = NULL,
  active_ratio_p = 0.02, 
  active_ratio_q = 0.005,
  residual_cor_mat = NULL,
  vec_rho_phenos = NULL,
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  50, 
  min_cor = 0, max_cor = 0.5,
  beta1 = 1,
  hm = 0.1,
  seed  = 2001
)

X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat


tol = 0.01
mu_t = 1
v_t = 4

devtools::load_all()
if(scheme == "Vanilla"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         tol = tol, 
                         CAVI = "full",
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
}else if(scheme == "RF"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         tol = tol, 
                         CAVI = "random",
                         start_partial = 50,
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
}else if (scheme == "AFE"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         tol = tol, 
                         CAVI = "adaptive",
                         start_partial = 50, 
                         epsilon_scheme = "ELBO",
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
}else if (scheme == "AFI"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         tol = tol, 
                         CAVI = "adaptive",
                         start_partial = 50, #Inf, 1000, 500, 100
                         epsilon_scheme = "iteration",
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
  
  
}else if(scheme =="AFIO"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         tol = tol, 
                         CAVI = "adaptive",
                         start_partial = 50, 
                         epsilon_scheme = "iteration",
                         eval_perform = T,
                         thinned_elbo_eval = T)
  })
}else if(scheme =="minimum"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         tol = tol, 
                         CAVI = "adaptive",
                         start_partial = 50, 
                         epsilon_scheme = "minimum",
                         min_epsilon = 0.1,
                         eval_perform = T,
                         thinned_elbo_eval = T)
  })
}
  
#Timing
perform_df = res_atlas$perform_df %>% 
  mutate(time_local = time_loop + time_tau+time_sig2_beta+time_m2_beta+time_Z+time_theta,
         time_global = time_theta + time_thetaPzeta + time_sigma) 

ROC.atlasqtl(res_atlas, pat_sim)

#visualize the partial updates
res_atlas$perform_df %>% 
  ggplot(aes(x = iter, y = subsample_size))+
  geom_point()
  
  