library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(ggplot2)
library(patchwork)
library(biomaRt)
library(stringr)
library(ggpubr)
library(gtools)

simulation_tool_path = "~/project/response_simulation"
source("~/project/response_simulation/data_simulation_tools/simulate_withRealGenome.R")
source("~/project/UKBB_run/utils/corr_plot_utils.R")
source("~/project/UKBB_run/utils/res_comb_utils.R")

save_path = "/rds/user/yl2021/hpc-work/hotspot_sim/simulation_res"
mybiobank = "/rds/user/yl2021/hpc-work/myukbb" 


task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

qt_ls = c(10, 50, 100, 250, 500, 1000) #set different number of active proteins
qt = qt_ls[task_id]


#################################################################################
#load data

#read snp_gene_info
snp_gene_info = fread("/rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chr19_filtered.bim")
colnames(snp_gene_info)[2] = "ID"
colnames(snp_gene_info)[4] = "POS"
snp_gene_info = snp_gene_info %>% dplyr::select(ID, POS) %>% 
  distinct(ID, .keep_all = TRUE) #remove duplicates

Y = readRDS("/rds/user/yl2021/hpc-work/myukbb/proteomics_clean/residuals_meanImp_discovery.rds") %>% as.data.table()
ID = Y$eid

# Y_real_corr = readRDS("/rds/user/yl2021/hpc-work/hotspot_sim/Y_real_corr.rds")
# protein_rank = readLines("/rds/user/yl2021/hpc-work/hotspot_sim/protein_rank.txt")

# geno = fread("/rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.raw") 
# X = geno[match(ID, geno$IID), .SD, .SDcols = c(2, 7:ncol(geno))]%>% 
#   setNames(map_chr(colnames(.), ~ strsplit(.x, "_")[[1]][1])) %>% 
#   dplyr::rename(eid = IID)
# X_real = as.matrix(X[,(ncol(X)-1000+1):ncol(X)])



################################################################################
# Data simulation


list = readRDS(file.path(save_path, paste0("dat/sim_list_", task_id, ".rds")))

X_sim = as.matrix(list$X)
Y_sim = as.matrix(list$Y)
pat_sim = list$pat
beta_sim = list$beta

p = 1000
q = 2919
n = nrow(X_sim)


################################################################################
# Run atlasQTL on simulated data
# the prior expectation and variance of the number of predictors associated with each response.
obj_atlasqtl = readRDS(file.path(save_path, paste0("atlasqtl_res/obj_atlasqtl_sim_", task_id, ".rds")))

# ------------------------------------------------------------------------------
# Infer FDR
# exact FDR
gam_vb = obj_atlasqtl$gam_vb
predicted <- gam_vb > 0.9

TP <- sum((predicted==1) & (pat_sim == 1))
FP <- sum(predicted & (pat_sim == 0))
FDR_exact <- FP / sum(predicted==1)


# loci-wise FDR
active_set_pred = which(apply(gam_vb, 2, function(col) any(col > 0.9)))
active_set_true = which(apply(pat_sim, 2, function(col) any(col == 1)))
TP_proteins = sum(active_set_pred %in% active_set_true)
FP_proteins = sum(!active_set_pred %in% active_set_true)

FDR_loci = FP_proteins / (FP_proteins + TP_proteins)

# AUROC, AUPRC
roc <- PRROC::roc.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUROC = roc$auc

pr <- PRROC::pr.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUPRC = pr$auc.integral


################################################################################
#run GWAS
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
  "--out /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim",
  "--glm"
)

system(plink_cmd_GWAS)

#################################################################################
#load GWAS res
gwas_res = read_PLINK_GWAS("/rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/", pattern = "chr19_sim\\.PHENO\\d+\\.glm\\.linear", protein_names=colnames(Y_sim))

# beta
beta_mat = as.data.frame(gwas_res[[2]])
beta_mat$ID = rownames(beta_mat)
beta_df = reshape2::melt(beta_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")

# pval
corr_mat = as.data.frame(gwas_res[[1]])
corr_mat$ID = rownames(corr_mat)

# joint 
gwas_df = reshape2::melt(corr_mat, id.vars = "ID", variable.name = "protein_name", value.name = "corr_metric")%>% 
  #Combine snp_gene_info
  left_join(snp_gene_info, by = "ID") %>% 
  left_join(beta_df, by = c("ID", "protein_name")) %>% 
  mutate(corr_metric = -log10(corr_metric))

#-------------------------------------------------------------------------------
# GWAS FDR

predicted <- gwas_res[[1]] < 0.05/(p*q)

TP <- sum((predicted==1) & (pat_sim == 1))
FP <- sum(predicted & (pat_sim == 0))
FDR_exact_GWAS <- FP / sum(predicted==1)


# loci-wise FDR
active_set_pred = which(apply(gwas_res[[1]], 2, function(col) any(col < 0.05/(p*q))))
active_set_true = which(apply(pat_sim, 2, function(col) any(col == 1)))
TP_proteins = sum(active_set_pred %in% active_set_true)
FP_proteins = sum(!active_set_pred %in% active_set_true)

FDR_loci_GWAS = FP_proteins / (FP_proteins + TP_proteins)

# # AUROC, AUPRC
# roc <- PRROC::roc.curve(scores.class0 = -log10(as.numeric(as.vector(as.matrix(gwas_res[[1]])))), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
# AUROC_GWAS = roc$auc
# 
# pr <- PRROC::pr.curve(scores.class0 = -log10(as.numeric(as.vector(as.matrix(gwas_res[[1]])))), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
# AUPRC_GWAS = pr$auc.integral



################################################################################
#make plots

reshape_atlasQTL_res = function(obj_atlasqtl){
  # beta
  beta_mat = as.data.frame(obj_atlasqtl$beta_vb)
  beta_mat$ID = rownames(beta_mat)
  beta_df = reshape2::melt(beta_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")
  
  # PPI
  corr_mat = as.data.frame(obj_atlasqtl$gam_vb)
  corr_mat$ID = rownames(corr_mat)
  
  # joint beta and ppi
  corr_df = reshape2::melt(corr_mat, id.vars = "ID", variable.name = "protein_name", value.name = "corr_metric")%>% 
    #Combine snp_gene_info
    left_join(snp_gene_info, by = "ID") %>% 
    left_join(beta_df, by = c("ID", "protein_name"))
  
  # theta
  theta_df = 
    data.frame(ID = names(obj_atlasqtl$theta_vb), theta = obj_atlasqtl$theta_vb) %>% 
    left_join(snp_gene_info, by = "ID")
  
  return(
    list(
      corr_df,theta_df
    )
  )
}

#-------------------------------------------------------------------------------
#manhattan plot of simulated v.s. inferred

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
res_simData = reshape_atlasQTL_res(obj_atlasqtl)
p1.2 = res_simData[[1]]%>% 
  group_by(POS, ID) %>% 
  summarise(n_assoc = sum(corr_metric > 0.9)) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  geom_point()+
  theme_bw()+
  ggtitle("Inferred # of associated proteins by SNP (atlasQTL)")


# GWAS inferred associations
p1.3 = gwas_df%>% 
  group_by(ID, POS) %>% 
  summarise(n_assoc = sum(corr_metric > -log10(0.05/(1000*3000)))) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  geom_point()+
  theme_bw()+
  ggtitle("Inferred # of associated proteins by SNP (GWAS)")

# #theta
# p1.3 = res_simData[[2]] %>% 
#   ggplot(aes(x = POS, y = theta))+
#   geom_point()+
#   theme_bw()+
#   ggtitle("Inferred theta")



p1 = wrap_plots(p1.1, p1.2, p1.3, ncol = 1)


#-------------------------------------------------------------------------------
# maximum PPI v.s. -log10pval
atlas_gwas_joint_df = res_simData[[1]] %>% left_join(gwas_df, by = c("ID", "POS","protein_name"))

p2.1  = atlas_gwas_joint_df %>% 
  group_by(protein_name) %>% 
  drop_na() %>% 
  summarise(max_ppi = max(corr_metric.x),
            max_logpval = max(corr_metric.y)) %>% 
  ggplot(aes(x =max_ppi, y = max_logpval))+
  xlab("maximum PPI by protein")+
  ylab(expression("maximum -log"[10]*"(p-value) by protein"))+
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = -log10(1.7e-8), linetype = "dashed", linewidth = 0.2)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2)+
  ggtitle("Maximum -log10pval v.s. PPI")+
  theme_bw()

# distributin of simulated BETA
p2.2 = data.frame(beta = as.vector(beta_sim)) %>% 
  filter(beta !=0) %>% 
  ggplot(aes(x = beta))+
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", alpha = 0.5) +
  ggtitle("Distribution of simulated BETA (non-zero)")+
  theme_bw()

p2 = wrap_plots(p2.2, p2.1, ncol = 1)
#-------------------------------------------------------------------------------
# Simulated BETA v.s. inferred
beta_sim_mat= as.data.frame(beta_sim)
beta_sim_mat$ID = rownames(beta_sim_mat)
beta_sim_df = reshape2::melt(beta_sim_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")

p3.1 = res_simData[[1]] %>% 
  left_join(beta_sim_df, by = c("ID", "protein_name")) %>% 
  group_by(protein_name) %>% 
  drop_na() %>% 
  summarise(max_BETA.x = max(abs(BETA.x)),
            max_BETA.y = max(abs(BETA.y)),
            max_ppi = max(corr_metric)) %>%
  ggplot(aes(x = abs(max_BETA.x), y = abs(max_BETA.y), color = max_ppi))+
  geom_point(alpha = 0.5)+
  xlab("BETA (atlasQTL)")+
  ylab("BETA (simulated)")+
  theme_bw()

p3.2 = gwas_df %>% 
  left_join(beta_sim_df, by = c("ID", "protein_name")) %>% 
  group_by(protein_name) %>% 
  drop_na() %>% 
  summarise(max_BETA.x = max(BETA.x),
            max_BETA.y = max(BETA.y),
            max_pval = max(corr_metric)) %>%
  ggplot(aes(x = abs(max_BETA.x), y = abs(max_BETA.y), color = max_pval))+
  geom_point(alpha = 0.5)+
  xlab("BETA (GWAS)")+
  ylab("BETA (simulated)")+
  theme_bw()

p3 = wrap_plots(p3.1, p3.2, ncol = 1)


# #max beta values, simulated v.s. inferred
# beta_sim_mat= as.data.frame(beta_sim)
# beta_sim_mat$ID = rownames(beta_sim_mat)
# beta_sim_df = reshape2::melt(beta_sim_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")
# 
# p2 = res_simData[[1]] %>% 
#   left_join(beta_sim_df, by = c("ID", "protein_name")) %>% 
#   group_by(protein_name) %>% 
#   drop_na() %>% 
#   summarise(max_BETA.x = max(abs(BETA.x)),
#             max_BETA.y = max(abs(BETA.y)),
#             max_ppi = max(corr_metric)) %>%
#   ggplot(aes(x = abs(max_BETA.x), y = abs(max_BETA.y), color = max_ppi))+
#   geom_point(alpha = 0.5)+
#   xlab("BETA (atlasQTL)")+
#   ylab("BETA (simulated)")+
#   theme_bw()+
#   ggtitle("inferred BETA v.s. simulated (maximum by protein)")


#save to subfolder plots
p = wrap_plots(p1, p2, p3, nrow = 1)
ggsave(filename = file.path(save_path, paste0("plots/plot_", task_id, ".pdf")), plot = p, device = "pdf", 
       width = 20, height = 8, dpi = 1200)




#-------------------------------------------------------------------------------
#save results
output =
  data.table(
    qt= qt,
    FDR_exact = FDR_exact,
    FDR_loci = FDR_loci,
    AUROC = AUROC,
    AUPRC = AUPRC,
    FDR_exact_GWAS = FDR_exact_GWAS,
    FDR_loci_GWAS = FDR_loci_GWAS
  )

write.csv(output, file.path(save_path, paste0("output_", task_id, '.csv')), row.names = FALSE)


