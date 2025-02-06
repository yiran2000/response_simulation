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
library(biomaRt)
library(stringr)

#change the file directory as in your case
mybiobank = "/rds/user/yl2021/hpc-work/myukbb"
geno_path = "/rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chunks"
run_path = "/rds/user/yl2021/hpc-work/UKB_run_results/chr19_2"
geno_files = gtools::mixedsort(list.files(geno_path,pattern = "\\.raw$", full.names = TRUE))
save_path = "/rds/user/yl2021/hpc-work/hotspot_sim"
simulation_tool_path = "~/project/response_simulation"

source("data_simulation_tools/simulate_withRealGenome.R")
source("data_simulation_tools/simulate_withRealGenome.R")
source("~/project/UKBB_run/utils/corr_plot_utils.R")
source("~/project/UKBB_run/utils/replication_rate_helpers.R")
source("~/project/UKBB_run/utils/permute_utils.R")
source("~/project/UKBB_run/utils/res_comb_utils.R")


#read snp_gene_info
snp_gene_info = fread("/rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chr19_filtered.bim")
colnames(snp_gene_info)[2] = "ID"
colnames(snp_gene_info)[4] = "POS"
snp_gene_info = snp_gene_info %>% dplyr::select(ID, POS) %>% 
  distinct(ID, .keep_all = TRUE) #remove duplicates


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




protein_rank = readLines(file.path(save_path, "protein_rank.txt"))
Y_real_corr = readRDS("/rds/user/yl2021/hpc-work/hotspot_sim/Y_real_corr.rds")


list = readRDS(file.path(save_path, "sim_list_3.rds"))

X_sim = as.matrix(list$X)
Y_sim = as.matrix(list$Y)
pat_sim = list$pat
beta_sim = list$beta

p = 1000
q = ncol(Y_real)
n = row(X_real)


obj_atlasqtl = readRDS("/rds/user/yl2021/hpc-work/hotspot_sim/obj_atlasqtl_sim_3.rds")


# manhattan plot
# simulated true associations
p1 = data.frame(ID = rownames(pat_sim), ProteinCount = rowSums(pat_sim)) %>% 
  left_join(snp_gene_info, by = "ID") %>% 
  ggplot(
    aes(x = POS, y = ProteinCount)
  ) +
  # scale_x_continuous(limits = gene_range)+
  geom_point()+
  theme_bw()

# atlasQTL inferred associations
res_simData = reshape_atlasQTL_res(obj_atlasqtl)
p2 = res_simData[[1]]%>% 
  group_by(POS, ID) %>% 
  summarise(n_assoc = sum(corr_metric > 0.9)) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  geom_point()+
  theme_bw()

p3 = res_simData[[2]] %>% 
  ggplot(aes(x = POS, y = theta))+
  geom_point()+
  theme_bw()

wrap_plots(p1, p2, p3, ncol = 1)


res_simData[[1]] %>% filter(ID == "rs3810353")

res_simData[[2]] %>% filter(POS>55681600, POS<55701600)%>% 
  ggplot(aes(x = POS, y = theta))+
  geom_vline(xintercept = 55691600)+
  geom_point()+
  theme_bw()

################################################################################
#summarise loci-wise FDR

# Infer FDR
# exact FDR
gam_vb = obj_atlasqtl$gam_vb
predicted <- gam_vb > 0.9

TP <- sum((predicted==1) & (pat_sim == 1))
FP <- sum(predicted & (pat_sim == 0))
FDR <- FP / sum(predicted==1)
FDR

# loci-wise FDR
active_set_pred = which(apply(gam_vb, 2, function(col) any(col > 0.9)))
active_set_true = which(apply(pat_sim, 2, function(col) any(col == 1)))
TP_proteins = sum(active_set_pred %in% active_set_true)
FP_proteins = sum(!active_set_pred %in% active_set_true)

FP_proteins / (FP_proteins + TP_proteins)

################################################################################
#a better joining strategy: consider distance and LD
distance = 1e5/10
ld_thres = 0.5

active_SNP_sim = rownames(pat_sim)[which(rowSums(pat_sim)>0)]
active_POS_sim = snp_gene_info %>% filter(ID %in% active_SNP_sim) %>% pull(POS) %>% sort()

diff(active_POS_sim)

#now brutally pull thresholds, use 1/2...
midpoints = (head(active_POS_sim, -1) + tail(active_POS_sim, -1)) / 2

#can we always find something nearby?




# #-------------------------------------------------------------------------------
# #Evaluate AUROC, AUPRC
roc <- PRROC::roc.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUROC = roc$auc

pr <- PRROC::pr.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUPRC = pr$auc.integral


#-------------------------------------------------------------------------------
#save X 
snp_list = unique(colnames(X_real))
writeLines(snp_list, file.path("/rds/user/yl2021/hpc-work/hotspot_sim/simulation_snplist.txt"))
plink_cmd <- paste(
  "plink2",
  "--bgen /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chr19_filtered.bgen ref-first",
  "--sample /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chr19_filtered.sample",
  "--chr 19",
  "--extract /rds/user/yl2021/hpc-work/hotspot_sim/simulation_snplist.txt",
  "--from-bp 54602362 --to-bp 55728080",
  "--rm-dup force-first",
  "--force-intersect",
  "--maf 0.01",
  "--keep /rds/user/yl2021/hpc-work/myukbb/filtered_samples/sample_qc_pass_id_discovery.txt",
  "--export bgen-1.3",
  "--out /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp"
)

system(plink_cmd)


#save Y
#re-shuffle Yto have the same order as X_real 
# row_mapping <- match(
#   apply(X_real, 1, paste, collapse = "_"),
#   apply(X_sim, 1, paste, collapse = "_")
# )
# 
# #Reorder Y_sim to match the order of X_real
# Y_sim_ordered <- Y_sim[row_mapping, ]



ID = fread("/rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/simulated_snp.fam")
pheno_sim = as.data.table(Y_sim) %>% 
  mutate(FID = Y$eid, IID = Y$eid) %>% 
  dplyr::select(FID, IID, everything())
pheno_sim_reordered = pheno_sim[match(ID$V1, pheno_sim$IID), ]

write.table(pheno_sim_reordered, file.path(mybiobank, "my_pheno/pheno_sim.fam"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)

#run GWAS
plink_cmd_GWAS <- paste(
  "plink2",
  "--bfile /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/simulated_snp",
  "--chr 19",
  "--extract /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/simulation_snplist.txt",
  "--keep /rds/user/yl2021/hpc-work/myukbb/filtered_samples/sample_qc_pass_id_discovery.txt",
  "--pheno /rds/user/yl2021/hpc-work/myukbb/my_pheno/pheno_sim.fam",
  "--out /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim",
  "--glm"
)

system(plink_cmd_GWAS)
# list.files("/rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/")

#-------------------------------------------------------------------------------
#LOAD GWAS results
library(gtools)
gwas_res = read_PLINK_GWAS("/rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/", pattern = "chr19_sim\\.PHENO\\d+\\.glm\\.linear", protein_names=colnames(Y_sim))
# gwas_df = process_corr_mat(corr_mat = as.data.frame(gwas_res[[1]]), beta_mat = as.data.frame(gwas_res[[2]]),
#                            protein_names = protein_names, snp_gene_info =snp_gene_info, protein_gene_info =protein_gene_info, 
#                            CHR=CHR, corr_thres=NULL)

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

#manhattan plot 
gwas_df%>% 
  group_by(ID, POS) %>% 
  summarise(n_assoc = sum(corr_metric > -log10(0.05/(1000*3000)))) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  geom_point()+
  theme_bw()




#compare maximum PPI and -log10pval 
atlas_gwas_joint_df = res_simData[[1]] %>% left_join(gwas_df, by = c("ID", "POS","protein_name"))

atlas_gwas_joint_df %>% 
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
  theme_bw()


#compare my inferred BETA with simulated values
beta_sim_mat= as.data.frame(beta_sim)
beta_sim_mat$ID = rownames(beta_sim_mat)
beta_sim_df = reshape2::melt(beta_sim_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")

p1 = res_simData[[1]] %>% 
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

p2 = gwas_df %>% 
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

p1+p2

