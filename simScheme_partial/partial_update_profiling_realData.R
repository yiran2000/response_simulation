library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(pryr)
library(pROC)
library(ggplot2)
library(patchwork)
library(egg)
library(cowplot)
library(splines)
library(atlasqtl)
library(tictoc)
library(purrr)
library(gtools)
library(peakRAM)
source("~/project/UKBB_run/utils/corr_plot_utils.R")
source("~/project/UKBB_run/utils/replication_rate_helpers.R")
source("~/project/UKBB_run/utils/permute_utils.R")
source("~/project/UKBB_run/utils/res_comb_utils.R")



mybiobank = "/rds/user/yl2021/hpc-work/myukbb"
geno_path = "/rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chunks"
run_path = "/rds/user/yl2021/hpc-work/UKB_run_results/chr19_selected_anneal"

geno_files = gtools::mixedsort(list.files(geno_path,pattern = "\\.raw$", full.names = TRUE))

task_id = if_else(!is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))), as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")), 1)

#-------------------------------------------------------------------------------
param_table = 
  expand.grid(
    chunk_id = c(3,8,17,28,30,37),
    method = c("full", "partial_eIter_adaptive_anneal"))

chunk_id = param_table[[task_id, "chunk_id"]]
method = param_table[[task_id, "method"]]


#-------------------------------------------------------------------------------
#

# Read in Y data
Y_discovery = readRDS(file.path(mybiobank, "/proteomics_clean/residuals_meanImp_discovery.rds")) %>% as.data.table()
ID_discovery = Y_discovery$eid


# Read in X data
geno_file_name = geno_files[chunk_id]
geno = fread(geno_file_name) 


X_discovery = geno[match(ID_discovery, geno$IID), .SD, .SDcols = c(2, 7:ncol(geno))]%>% 
  setNames(map_chr(colnames(.), ~ strsplit(.x, "_")[[1]][1])) %>% 
  dplyr::rename(eid = IID)

#Prepare 
q = ncol(Y_discovery)-2
p = ncol(X_discovery)-1

#Define mu_t, v_t
# mu_t = 0.008 * q/p
# v_t = 0.06 * (q/p)^2

# mu_t = 2e-5* p
# std = 2e-4*p
# v_t = max(mu_t, std^2)

mu_t = 1
v_t = 4

start = extract_start_position(geno_file_name)
end = extract_end_position(geno_file_name)
  
#-------------------------------------------------------------------------------
#Run atlasQTL

if(method == "full"){
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_discovery[,3:ncol(Y_discovery)]), X = as.matrix(X_discovery[,2:ncol(X_discovery)]),
                         p0 = c(mu_t, v_t),
                         user_seed = 2000, maxit=3000,
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
}else{
  time = peakRAM({
    res_atlas = atlasqtl(Y = as.matrix(Y_discovery[,3:ncol(Y_discovery)]), X = as.matrix(X_discovery[,2:ncol(X_discovery)]),
                         p0 = c(mu_t, v_t),
                         user_seed = 2000, maxit=3000,
                         batch = "y",
                         tol = 0.01, 
                         start_partial = 10000, #Inf, 1000, 500, 100
                         max_partial = Inf,
                         end_full = F,
                         subsample_scheme = "geometric",
                         min_epsilon = 0,
                         geom_alpha = 0.95,
                         eval_perform = T,
                         thinned_elbo_eval = T)
  })
}


#-------------------------------------------------------------------------------
#save results
output = 
  data.table(
    chunk = chunk_id,
    method = method,
    start = start,
    end = end,
    p = p,
    n_iter = res_atlas$it,
    runtime_min = time$Elapsed_Time_sec/60,
    memory_peak = time$Peak_RAM_Used_MiB,
    memory_total = time$Total_RAM_Used_MiB,
    ELBO = res_atlas$lb_opt
  )

write.csv(output, file.path(run_path, paste0("output_", task_id, '.csv')), row.names = FALSE)
saveRDS(res_atlas, file.path(run_path, paste0("obj_atlasqtl_", task_id, '.rds')))

