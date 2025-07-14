library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(tictoc)
library(PRROC)

library(atlasqtl)
library(echoseq)

library(peakRAM)

################################################################################
#load functions

HPC = T
task_id = if_else(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))), as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")), 1)


if(HPC){
  simulation_tool_path = "~/project/response_simulation"
  list = readRDS(file.path("/rds/user/yl2021/hpc-work/hotspot_sim", "sim_list_3.rds"))
  X_real = as.matrix(list$X)
  save_path = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme_3"
  
}else{
  simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
  #Load in Real SNP data
  list = readRDS(file.path(simulation_tool_path, "data/sim_list_3.rds"))
  X_real = as.matrix(list$X)
}

source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source(file.path(simulation_tool_path, "demo_tools/res_demo_helpers.R"))


rm(list)
gc()

#First compare full & partial_eELBO, partial_eIter, with starting point = 500 (5)
#And then in a specific setting, compare different parameter choices

################################################################################
#Define parameters


full_params <- expand.grid(
  rep_id = 1:25,
  scheme = "full"
) %>% mutate(start_partial = Inf, min_epsilon =0)  

partial_params <- expand.grid(
  rep_id = 1:25,
  start_partial = c(1000, 500, 100),
  min_epsilon = c(0, 0.01, 0.1),
  scheme = c("partial_eELBO", "partial_eIter_0.99", "partial_eIter_0.95", "partial_eIter_0.9","partial_eIter_0")
)


param_table <- bind_rows(full_params, partial_params)

hm = 0.15
active_ratio_p = 0.003
active_ratio_q = 0.005
q = 1000
p = 3000
mu_t = 1
v_t = 4
start_partial = param_table[[task_id, "start_partial"]]
seed = param_table[[task_id, "rep_id"]]
scheme = param_table[[task_id, "scheme"]]
min_epsilon = param_table[[task_id, "min_epsilon"]]
################################################################################
#simulate one single hotspot
list_sim = simulate_withRealGenome(
  X = X_real[1:p,], 
  q = q,
  protein_ls = paste0("protein_", 1:q),
  active_protein = NULL,
  active_ratio_p = active_ratio_p, 
  active_ratio_q = active_ratio_q,
  residual_cor_mat = NULL,
  vec_rho_phenos = NULL,
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  50, 
  min_cor = 0, max_cor = 0.5,
  hm = hm,
  seed  = seed
)


X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat

################################################################################
#run atlasQTL

if(scheme == "full"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = 1,
                         end_full = F,
                         epsilon_scheme = "geometric",
                         min_epsilon = min_epsilon,
                         geom_alpha = 0.97,
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
}else if(scheme == "partial_eELBO"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         epsilon_scheme = "logistic",
                         min_epsilon = min_epsilon,
                         geom_alpha = 0.97,
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
  
}else if (scheme == "partial_eIter_0.99"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         epsilon_scheme = "geometric",
                         min_epsilon = min_epsilon,
                         geom_alpha = 0.99,
                         eval_perform = T,
                         thinned_elbo_eval = T)
  })
}else if (scheme == "partial_eIter_0.95"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         epsilon_scheme = "geometric",
                         min_epsilon = min_epsilon,
                         geom_alpha = 0.95,
                         eval_perform = T,
                         thinned_elbo_eval = T)
  })
}else if (scheme == "partial_eIter_0.9"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         epsilon_scheme = "geometric",
                         min_epsilon = min_epsilon,
                         geom_alpha = 0.9,
                         eval_perform = T,
                         thinned_elbo_eval = T)
    
  })
}else if (scheme == "partial_eIter_0"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         epsilon_scheme = "geometric",
                         min_epsilon = min_epsilon,
                         geom_alpha = 0,
                         eval_perform = T,
                         thinned_elbo_eval = T)
    
  })
}
# res_atlas$perform_df %>% ggplot(aes(x = iter, y = e, color =partial))+geom_point()

################################################################################
# iterations
summary_df <- data.frame(
  task_id = task_id, 
  rep_id = param_table[[task_id, "rep_id"]],
  method = scheme,
  start_partial = start_partial,
  min_epsilon = min_epsilon,
  iterations = res_atlas$it, 
  iterations_full = sum(res_atlas$perform_df$partial==T),
  runtime = time$Elapsed_Time_sec,
  memory_peak = time$Peak_RAM_Used_MiB,
  memory_total = time$Total_RAM_Used_MiB,
  ELBO = res_atlas$lb_opt, 
  AUROC = ROC.atlasqtl(res_atlas, pat_sim)$AUROC,
  AUPRC = ROC.atlasqtl(res_atlas, pat_sim)$AUPRC
)

write.table(summary_df,file.path(save_path, paste0("summary/summary_", task_id, ".csv")))
write.table(res_atlas$theta_vb, file.path(save_path, paste0("theta/theta_", task_id, ".csv")), col.names = F)
write.table(res_atlas$zeta_vb,file.path(save_path, paste0("zeta/zeta_", task_id, ".csv")), col.names = F)
write.table(res_atlas$gam_vb,file.path(save_path, paste0("gam/gam_", task_id, ".csv")))



