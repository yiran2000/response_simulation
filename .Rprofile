# rm()
# gc()
# 
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(data.table)
# library(ggplot2)
# library(atlasqtl)
# library(tictoc)
# library(echoseq)
# library(PRROC)
# library(patchwork)
# library(stringr)
# library(gtools)
# 
# 
# simulation_tool_path = "~/project/response_simulation"
# data_path = "/rds/user/yl2021/hpc-work/hotspot_sim"
# save_path = "/rds/user/yl2021/hpc-work/hotspot_sim/simulation_res"
# mybiobank = "/rds/user/yl2021/hpc-work/myukbb" 
# 
# # simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
# # data_path = "data"
# 
# 
# source("~/project/response_simulation/data_simulation_tools/simulate_withRealGenome.R")
# source("~/project/response_simulation/reshape_atlasQTL_res.R")
# 
# source("~/project/UKBB_run/utils/corr_plot_utils.R")
# source("~/project/UKBB_run/utils/res_comb_utils.R")
# 
# 
# ################################################################################
# #read in pre-save data
# #X
# list = readRDS(file.path(data_path, "sim_list_3.rds"))
# X_real = as.matrix(list$X)
# 
# rm(list)
# 
# Y = readRDS("/rds/user/yl2021/hpc-work/myukbb/proteomics_clean/residuals_meanImp_discovery.rds") %>% as.data.table()
# Y_real_corr = readRDS(file.path(data_path,"Y_real_corr.rds"))
# protein_rank = readLines(file.path(data_path,"protein_rank.txt"))
# snp_gene_info = fread(file.path(data_path,"snp_gene_info.csv"))

