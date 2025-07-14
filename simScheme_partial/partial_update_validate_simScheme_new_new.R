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
  save_path = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/Sim1_new_new"
  
  # if(task_id == 1){
  #   system("find /rds/user/yl2021/hpc-work/slurmdir/ -type f -delete")
  #   system("find /rds/user/yl2021/hpc-work/partial_update_validate_simScheme/Sim1_new_new/ -type f -delete")
  # }
  
  
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
    data_id = 1:60,
    active_ratio_q = c(0.005, 0.2, 0.5), #it is ok to have 0.05 but not really neccessary
    hm = c(0.05, 0.15, 0.3))
param_table$task_id = 1:nrow(param_table)


param_table = 
  expand.grid(
    data_id = 1:90,
    active_ratio_q = c(0.005, 0.2, 0.5), #it is ok to have 0.05 but not really neccessary
    hm = c(0.05, 0.3, 0.15))
param_table$task_id = 1:nrow(param_table)

param_table %>% filter(hm == 0.15, active_ratio_q== 0.5)


# param_table %>% filter(active_ratio_q == 0.5, hm == 0.3)

q = 3000
p = 1000
active_ratio_p = 0.01
tol = 0.01

# active_ratio_q = 0.005
# hm = 0.3
# data_id = 2

data_id = param_table[[task_id, "data_id"]]
scheme = param_table[[task_id, "scheme"]]
active_ratio_q = param_table[[task_id, "active_ratio_q"]]
hm = param_table[[task_id, "hm"]]

# if(param_table[[task_id, "case"]] =="sparse"){
#   active_ratio_p = 0.002
#   active_ratio_q = 0.002 
# }else{
#   active_ratio_p = 0.005
#   active_ratio_q = 0.2 
# }

################################################################################
#simulate one single hotspot
list_sim = simulate_withRealGenome(
  X = X_real[,1:p], 
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

for(scheme in c("Vanilla", "RF","AFE", "AFI","AFIO")){
  if(scheme == "Vanilla"){
    time = peakRAM({
      res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                           p0 = c(mu_t, v_t),
                           user_seed = 25, maxit=3000,
                           batch = "y",
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
                           batch = "y",
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
                           batch = "y",
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
                           batch = "y",
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
                           batch = "y",
                           tol = tol, 
                           CAVI = "adaptive",
                           start_partial = 50, 
                           epsilon_scheme = "iteration",
                           eval_perform = T,
                           thinned_elbo_eval = T)
    })
  }
  
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
  
  #Timing
  perform_df = res_atlas$perform_df %>% 
    mutate(time_local = time_loop + time_tau+time_sig2_beta+time_m2_beta+time_Z+time_theta,
           time_global = time_theta + time_thetaPzeta + time_sigma) 
  
  #summary table
  summary_df <- data.frame(
    task_id = task_id,
    data_id = data_id,
    method = scheme,
    hm = hm,
    tol = tol,
    active_ratio_q = active_ratio_q,
    iterations = res_atlas$it, 
    runtime = time$Elapsed_Time_sec,
    time_init = res_atlas$time_init,
    time_loop = sum(perform_df$time_loop),
    time_local = sum(perform_df$time_local),
    time_global = sum(perform_df$time_global),
    time_total = sum(perform_df$time_total),
    time_elbo = sum(perform_df$time_ELBO),
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
  
  #save
  write.table(summary_df,file.path(save_path, paste0("summary/summary_", task_id, "_", scheme, ".csv")))
  write.table(perform_df,file.path(save_path, paste0("perform/perform_", task_id, "_", scheme, ".csv")))
  
  #clear up disk
  rm(res_atlas)
  rm(perform_df)
  rm(summary_df)
  gc()
}



