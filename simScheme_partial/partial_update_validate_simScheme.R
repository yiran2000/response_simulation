library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(tictoc)
library(PRROC)

library(atlasqtl)
library(echoseq)

################################################################################
#load functions

HPC = F

if(HPC){
  simulation_tool_path = "~/project/response_simulation"
}else{
  simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
}

source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source(file.path(simulation_tool_path, "demo_tools/res_demo_helpers.R"))


#Load in Real SNP data
list = readRDS(file.path(simulation_tool_path, "data/sim_list_3.rds"))
X_real = as.matrix(list$X)


################################################################################
#Define parameters
param_table = 
  expand.grid(rep_id = 1:20,
              scheme = c("full", "partial", "partial_full"),
              case = c("sparse","dense"),
              q = c(1000, 3000),
              hm = c(0.1,0.3))

seed = param_table[[task_id, "rep_id"]]*100
scheme = param_table[[task_id, "scheme"]]
q = param_table[[task_id, "q"]]
hm = param_table[[task_id, "hm"]]
if(param_table[[task_id, "case"]] =="sparse"){
  active_ratio_p = 0.002
  active_ratio_q = 0.002 
}else{
  active_ratio_p = 0.005
  active_ratio_q = 0.2 
}

################################################################################
#simulate one single hotspot
list_sim = simulate_withRealGenome(
  X = X_real, 
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
mu_t = 1
v_t = 4

if(scheme == "full"){
  time = system.time({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                              p0 = c(mu_t, v_t),
                              user_seed = 1, maxit=1000,
                              batch = "y",
                              tol = 0.01, 
                              max_partial = 1,
                              end_full = T,
                              epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                              eval_perform = T,
                              thinned_elbo_eval = T)
  })
}else if(scheme == "partial"){
  time = system.time({
    res_atlas_partial = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                                 p0 = c(mu_t, v_t),
                                 user_seed = 1, maxit=1000,
                                 batch = "y",
                                 tol = 0.01, 
                                 max_partial = Inf,
                                 end_full = F,
                                 epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                                 eval_perform = T,
                                 thinned_elbo_eval = T)
  })
  
}else{
  time = system.time({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                                      p0 = c(mu_t, v_t),
                                      user_seed = 1, maxit=1000,
                                      batch = "y",
                                      tol = 0.01, 
                                      max_partial = Inf,
                                      end_full = T,
                                      epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                                      eval_perform = T,
                                      thinned_elbo_eval = T)
  })
}


################################################################################
# iterations
summary_df <- data.frame(
  Repeat = param_table[[task_id, "rep_id"]],
  Method = scheme,
  Iterations = res_atlas$it, 
  runtime = time[[3]],
  ELBO = res_atlas$lb_opt, 
  AUROC = ROC.atlasqtl(res_atlas, pat_sim)$AUROC,
  AUPRC = ROC.atlasqtl(res_atlas, pat_sim)$AUPRC
)

writeLines(res_atlas$theta_vb, file.path(save_path, paste0("theta/theta_", task_id, ".csv")))
writeLines(res_atlas$zeta_vb,file.path(save_path, paste0("zeta/zeta_", task_id, ".csv")))
write.table(res_atlas$gam_vb,file.path(save_path, paste0("gam/gam_", task_id, ".csv")))
write.table(summary_df,file.path(save_path, paste0("summary/summary_", task_id, ".csv")))

