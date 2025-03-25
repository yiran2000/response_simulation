rm()
gc()

library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(patchwork)
library(stringr)


simulation_tool_path = "~/project/response_simulation"
data_path = "/rds/user/yl2021/hpc-work/hotspot_sim"
save_path = "/rds/user/yl2021/hpc-work/hotspot_sim/simulation_res"

# simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
# data_path = "data"


source("~/project/response_simulation/data_simulation_tools/simulate_withRealGenome.R")
source("~/project/response_simulation/reshape_atlasQTL_res.R")




################################################################################
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = 1))

# rep = 1:50
qt_ls = c(10, 50, 500, 1000) #set different number of active proteins
m_ls = c(0.01, 0.005, 0.001)

param_table = expand.grid(rep_id = 1:100, qt = qt_ls, m = m_ls)
qt = param_table[task_id,"qt"]
m = param_table[task_id,"m"]
rep_id = param_table[task_id,"rep_id"]

################################################################################
#read in pre-save data
#X
list = readRDS(file.path(data_path, "sim_list_3.rds"))
X_real = as.matrix(list$X)

rm(list)

Y_real_corr = readRDS(file.path(data_path,"Y_real_corr.rds"))
protein_rank = readLines(file.path(data_path,"protein_rank.txt"))
snp_gene_info = fread(file.path(data_path,"snp_gene_info.csv"))


################################################################################
#simulate one single hotspot
p = 200
L <- chol(Y_real_corr)

list_sim = simulate_withRealGenome(
    X = X_real[,801:1000], 
    q = 2919,
    protein_ls = colnames(Y_real_corr),
    active_SNP = NULL,
    active_protein = protein_rank[1:qt],
    active_ratio_p = 1/p, 
    active_ratio_q = NULL,
    residual_cor_mat = NULL,
    L = L,
    missing_ratio = 0, 
    max_tot_pve = NULL, 
    m = m,
    sh2 = 5,
    seed  = task_id
  )

#equvicorrelated
list_sim = simulate_withRealGenome(
  X = X_real[,801:1000], 
  q = 400,
  protein_ls = paste0("protein_", 1:400),
  active_protein = paste0("protein_", 1:150),
  active_ratio_p = 1/200, 
  active_ratio_q = NULL,
  residual_cor_mat = NULL,
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  200, 
  min_cor = 0.9, max_cor = 0.95,
  m = 0.001,
  sh2 = 5,
  seed  = task_id
)


# saveRDS(list_sim, "/rds/user/yl2021/hpc-work/hotspot_sim/list_sim_5.rds")
# list_sim = readRDS("/rds/user/yl2021/hpc-work/hotspot_sim/list_sim_4.rds")

X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat

# Y_sim_corr = cor(Y_sim)
hist(as.vector(beta_sim)[as.vector(beta_sim)!=0])

max_beta = max(as.vector(beta_sim))

################################################################################
#run atlasQTL

# 
# prior expectation and variance of the number of predictors associated with each response.
# num_assoc = colSums(pat_sim)
# mu_t = max(0.05, mean(num_assoc))
# v_t = max(1, var(num_assoc)*5)


mu_t = 1
v_t = 4

devtools::load_all()
obj_atlasqtl =  atlasqtl(Y = Y_sim, X = X_sim,
                         p0 = c(mu_t, v_t),
                         user_seed = 1, 
                         maxit = 3000,
                         batch = "y",
                         tol_loose = 1,
                         tol_tight = 0.1,
                         burn_in = 10,
                         maxit_full = 5,
                         maxit_subsample = 10,
                         n_partial_update = 500,
                         iter_ladder = c(5, 10, 15, 20, 25, 30, 40, 60, 80, 100),
                         e_ladder = c(0.9, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05),
                         eval_perform = F,
                         thinned_elbo_eval = T) 

################################################################################
#-------------------------------------------------------------------------------
gam_vb = obj_atlasqtl$gam_vb
res_simData = reshape_atlasQTL_res(obj_atlasqtl, snp_gene_info)
# rm(obj_atlasqtl)

hist(as.vector(gam_vb)[as.vector(gam_vb)>0.1], breaks = 100)
#-------------------------------------------------------------------------------
#manhattan plot of simulated v.s. inferred

if(rep_id == 1){
  # simulated true associations
  p1.1 = data.frame(ID = rownames(pat_sim), ProteinCount = rowSums(pat_sim)) %>% 
    left_join(snp_gene_info, by = "ID") %>% 
    ggplot(
      aes(x = POS, y = ProteinCount)
    ) +
    geom_point()+
    theme_bw()+
    ggtitle("Simulated # of associated proteins by SNP")
  
  # atlasQTL inferred associations
  p1.2 = res_simData[[1]]%>% 
    group_by(POS, ID) %>% 
    summarise(n_assoc = sum(corr_metric > 0.9)) %>% 
    ggplot(aes(x = POS, y = n_assoc))+
    geom_point()+
    theme_bw()+
    ggtitle("Inferred # of associated proteins by SNP (atlasQTL)")
  
  p1.3 = res_simData[[2]] %>% 
    ggplot(aes(x = POS, y = theta))+
    geom_point()+
    theme_bw()+
    ggtitle("Inferred theta (atlasQTL)")
  
  p1 = wrap_plots(p1.1, p1.2, p1.3, ncol = 1)
  ggsave(filename = file.path(save_path, paste0("plots/plot_", task_id, ".pdf")), plot = p1, device = "pdf", 
         width = 5, height = 9, dpi = 1200)
}

which(colSums(pat_sim) > 0)


#-------------------------------------------------------------------------------
# exact FDR
predicted <- gam_vb > 0.9

TP_exact <- sum((predicted==1) & (pat_sim == 1))
FP_exact <- sum(predicted & (pat_sim == 0))

TN_exact = sum((predicted==0) & (pat_sim == 0))
FN_exact = sum((predicted==0) & (pat_sim == 1))

FPR_exact <- FP_exact / sum(pat_sim == 0)
FPR_exact

FDR_exact <- FP_exact / sum(predicted==1)
FDR_exact

power_exact = TP_exact/sum(pat_sim == 1)
power_exact

#-------------------------------------------------------------------------------
# loci-wise FDR
active_set_pred = which(apply(gam_vb, 2, function(col) any(col > 0.9)))
active_set_true = which(colSums(pat_sim)>0)

TP_protein = sum(active_set_pred %in% active_set_true)
FP_protein = sum(!active_set_pred %in% active_set_true)

FN_protein = sum(!active_set_true %in% active_set_pred)
TN_protein = ncol(gam_vb) - (TP_protein + FP_protein + FN_protein)

FPR_protein = FP_protein / (FP_protein + TN_protein)
FDR_protein = FP_protein / (TP_protein + FP_protein)
power_protein = TP_protein / (TP_protein + FN_protein)


roc <- PRROC::roc.curve(scores.class0 = as.vector(gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUROC = roc$auc

pr <- PRROC::pr.curve(scores.class0 = as.vector(gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUPRC = pr$auc.integral

################################################################################
#save results


output = data.table(
  task_id = task_id,
  m = m,
  qt = qt,
  rep_id = rep_id,
  max_beta = max_beta,
  AUROC = AUROC,
  AUPRC = AUPRC,
  TP_exact = TP_exact,
  FP_exact = FP_exact,
  TN_exact = TN_exact,
  FN_exact = FN_exact,
  FDR_exact = FDR_exact,
  FPR_exact = FPR_exact,
  power_exact = power_exact,
  TP_protein = TP_protein,
  FP_protein = FP_protein,
  TN_protein = TN_protein,
  FN_protein = FN_protein,
  FDR_protein = FDR_protein,
  FPR_protein = FPR_protein,
  power_protein = power_protein
)

write.csv(output, file.path(save_path, paste0("output_", task_id, '.csv')), row.names = FALSE)


#################################################################################
#save current snplist for GWAS 
snplist = colnames(X_sim)
write.table(snplist, file = file.path(data_path, "snplist.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#save simulated Y
ID = fread("/rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample")
pheno_sim = as.data.table(Y_sim) %>% 
  mutate(FID = Y$eid, IID = Y$eid) %>% 
  dplyr::select(FID, IID, everything())
pheno_sim_reordered = pheno_sim[match(ID$ID_1, pheno_sim$IID), ]

write.table(pheno_sim_reordered, file.path(mybiobank, "my_pheno/pheno_sim.fam"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)

#run GWAS
plink_cmd_GWAS <- paste(
  "plink2",
  "--bgen /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.bgen ref-first",
  "--sample /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample",
  "--chr 19",
  "--pheno /rds/user/yl2021/hpc-work/myukbb/my_pheno/pheno_sim.fam",
  "--extract /rds/user/yl2021/hpc-work/hotspot_sim/snplist.txt",
  "--out /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim",
  "--glm"
)

system(plink_cmd_GWAS)
