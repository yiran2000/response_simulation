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
  save_path = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/Sim3"
  
  if(task_id == 1){
    if(dir.exists("/rds/user/yl2021/hpc-work/slurmdir")) unlink("/rds/user/yl2021/hpc-work/slurmdir", recursive = TRUE)
    dir.create("/rds/user/yl2021/hpc-work/slurmdir")
    system("find /rds/user/yl2021/hpc-work/partial_update_validate_simScheme/Sim3/ -type f -delete")
  }
  

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
# Final simulation
# fix active ratio p
# 3 active ratio q
# 3 hm levels 
# 3 n
# 3 p
# 3 q
# full vs best partial 
# rep 15 times
# 960 jobs in total

#TODO
# also find a way to report the "running memory"
# test out the mu_t, v_t stuff
param_table = 
  expand.grid(
    data_id = 1:50,
    method = c("Vanilla", "AFIO"),
    active_ratio_q = c(0.001, 0.1, 0.02, 0.05, 0.1, 0.25, 0.5))



q = 3000
active_ratio_p = 0.005 # 2/200
hm = 0.1
start_partial = 1000

data_id = param_table[[task_id, "data_id"]]
method = param_table[[task_id, "method"]]
p = param_table[[task_id, "p"]]
active_ratio_q = param_table[[task_id, "active_ratio_q"]]


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
  seed  = data_id
)


X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat



################################################################################
#run atlasQTL
mu_t = 1
v_t = 4

if(method == "Vanilla"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = Inf, #Inf, 1000, 500, 100
                         max_partial = 1,
                         end_full = F, 
                         subsample_scheme = "geometric",
                         min_epsilon = 0,
                         geom_alpha = 0.97,
                         eval_perform = T,
                         thinned_elbo_eval = F)
  })
}else if (method == "AFIO"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 25, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = start_partial, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         subsample_scheme = "geometric",
                         min_epsilon = 0,
                         geom_alpha = 0.95,
                         eval_perform = T,
                         thinned_elbo_eval = T)
  })
}


# res_atlas$perform_df %>% ggplot(aes(x = iter, y = ELBO, color =partial))+geom_point()
# plot(res_atlas)
# plot(x = 1:200,y = rowSums(pat_sim ==1))
# res_atlas$zeta_ls

# ggplot(posterior_fittingCurve, aes(x = iter, y = value, color = associated)) +
#   geom_point(size = 0.5) +
#   scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
#   theme_minimal()
################################################################################
threshold <- 0.5

scores <- as.vector(res_atlas$gam_vb)  # predicted scores
labels <- as.numeric(as.vector(as.matrix(pat_sim)))  # ground truth (0 or 1)
predicted <- ifelse(scores >= threshold, 1, 0)
# caret::confusionMatrix(predicted, labels)
# Confusion matrix components
TP <- sum(predicted == labels & labels == 1)
TN <- sum(predicted == labels & labels == 0)
FN <- sum(predicted != labels & labels == 1)
FP <- sum(predicted != labels & labels == 0)

summary_df <- data.frame(
  task_id = task_id, 
  data_id = data_id,
  method = method,
  hm = hm,
  q = q,
  active_ratio_q = active_ratio_q,
  iterations = res_atlas$it, 
  iterations_full = sum(res_atlas$perform_df$partial==T),
  runtime = time$Elapsed_Time_sec,
  memory_peak = time$Peak_RAM_Used_MiB,
  memory_total = time$Total_RAM_Used_MiB,
  ELBO = res_atlas$lb_opt, 
  AUROC = ROC.atlasqtl(res_atlas, pat_sim)$AUROC,
  AUPRC = ROC.atlasqtl(res_atlas, pat_sim)$AUPRC,
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN
)

write.table(summary_df,file.path(save_path, paste0("summary/summary_", task_id, ".csv")))

# write.table(res_atlas$theta_vb, file.path(save_path, paste0("theta/theta_", task_id, ".csv")), col.names = F)
# write.table(res_atlas$zeta_vb,file.path(save_path, paste0("zeta/zeta_", task_id, ".csv")), col.names = F)
# write.table(res_atlas$gam_vb,file.path(save_path, paste0("gam/gam_", task_id, ".csv")))
# write.table(res_atlas$perform_df,file.path(save_path, paste0("perform/perform_", task_id, ".csv")))
# 
# #fitting curves
# zeta_fittingCurve = reshape_fittingCurve(res_atlas$zeta_ls, pat_sim)
# posterior_fittingCurve = reshape_fittingCurve(res_atlas$r_vc_ls, pat_sim)
# write.table(zeta_fittingCurve, file.path(save_path, paste0("zeta_fittingCurve/zeta_fittingCurve_", task_id, ".csv")))
# write.table(zeta_fittingCurve, file.path(save_path, paste0("posterior_fittingCurve/posterior_fittingCurve_", task_id, ".csv")))
