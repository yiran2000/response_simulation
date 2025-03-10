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

#filter and save snp_gene_info for those in this chunk
snp_gene_info = snp_gene_info %>% filter(ID %in% colnames(X))
write.csv(snp_gene_info, "/rds/user/yl2021/hpc-work/hotspot_sim/snp_gene_info.csv", row.names = FALSE)


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
# load Y
Y = readRDS(file.path(mybiobank, "/proteomics_clean/residuals_meanImp_discovery.rds")) %>% as.data.table()
ID = Y$eid

#-------------------------------------------------------------------------------
#Load X from GP6 loci
#PLINK chunk 37 remove duplicates and maf < 0.01

geno = fread("/rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.raw") 
X = geno[match(ID, geno$IID), .SD, .SDcols = c(2, 7:ncol(geno))]%>% 
  setNames(map_chr(colnames(.), ~ strsplit(.x, "_")[[1]][1])) %>% 
  dplyr::rename(eid = IID)

# #recode A into .raw file
# plink_cmd <- paste(
#   "plink2",
#   "--bgen /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chr19_filtered.bgen ref-first",
#   "--sample /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/chr19_filtered.sample",
#   "--chr 19",
#   "--extract /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/plink2.prune.in",
#   "--from-bp 54602362 --to-bp 55728080",
#   "--rm-dup force-first",
#   "--force-intersect",
#   "--maf 0.01",
#   "--keep /rds/user/yl2021/hpc-work/myukbb/filtered_samples/sample_qc_pass_id_discovery.txt",
#   "--export A",
#   "--out /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp"
# )
# #"--extract /rds/user/yl2021/hpc-work/myukbb/UKB_run_prepare/chr19/simulation_snplist.txt"
# system(plink_cmd)





#just take the range of one gene, and then simulate only 1 actual causal SNP associ. with 50 proteins, 
#add normal level of correlation/correlation from actual data
# mart = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl", GRCh=37,host = "http://www.ensembl.org") #UKBB uses ,GRCh=37
# GP6_range = getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"),
#   filters = c("chromosome_name", "external_gene_name"),
#   values = list("19", "GP6"),
#   mart = mart
# )
# gene_range  = c(GP6_range$start_position, GP6_range$end_position)

# gene_range = c(55525073,55549632)
# 
# snp_gene_info = fread(file.path(mybiobank, "UKB_run_prepare/chr19/chr19_filtered.bim"))
# colnames(snp_gene_info)[2] = "ID"
# colnames(snp_gene_info)[4] = "POS"
# snp_gene_info = snp_gene_info %>% dplyr::select(ID, POS) %>% 
#   distinct(ID, .keep_all = TRUE) #remove duplicates
# 
# #filter out the SNPs in the GP6 loci
# GP6_SNP = snp_gene_info %>% filter(POS>gene_range[1], POS< gene_range[2]) %>% pull(ID) %>% unique()
# X_GP6 = X[,..GP6_SNP]


#-------------------------------------------------------------------------------
# Apply atlasQTL on the actual chunk of Y and obtain the residuals
Y_real = as.matrix(Y[,3:ncol(Y)])
X_real = as.matrix(X[,(ncol(X)-1000+1):ncol(X)])



p = ncol(X_real)
q = ncol(Y_real)
mu_t = 2e-5* p
std = 2e-4*p
v_t = max(mu_t, std^2)


obj_atlasqtl_realData =  atlasqtl(Y = Y_real, X = X_real,
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=3000,
                         batch = "y",
                         tol_loose = 1,
                         tol_tight = 0.1,
                         burn_in = 10,
                         maxit_full = 5000,
                         maxit_subsample = 10,
                         n_partial_update = 500,
                         iter_ladder = c(5, 10, 15, 20, 25, 30, 40, 60, 80, 100),
                         e_ladder = c(0.9, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05),
                         eval_perform = F,
                         thinned_elbo_eval = T)

#plot results
res_realData = reshape_atlasQTL_res(obj_atlasqtl_realData)
#gene_range = range(res_realData[[1]]$POS)
p0 = res_realData[[1]]%>% 
  group_by(POS, ID) %>% 
  summarise(n_assoc = sum(corr_metric > 0.9)) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  scale_x_continuous(limits = gene_range)+
  geom_point()+
  theme_bw()

res_realData[[2]] %>% 
  ggplot(aes(x = POS, y = theta))+
  geom_point()+
  theme_bw()


#extract top 500 proteins by PPI
active_protein = res_realData[[1]] %>% group_by(protein_name) %>% summarise(max_ppi = max(corr_metric)) %>% arrange(desc(max_ppi)) %>% pull(protein_name) %>% head(500) %>% as.character()
protein_rank = res_realData[[1]] %>% group_by(protein_name) %>% summarise(max_ppi = max(corr_metric)) %>% arrange(desc(max_ppi)) %>% pull(protein_name) %>% as.character()
writeLines(protein_rank, file.path(save_path, "protein_rank.txt"))

# #distribution of beta
# res_summ_protein = res_realData[[1]]%>% 
#   group_by(protein_name) %>% 
#   summarise(max_BETA = max(BETA),
#             max_PPI = max(corr_metric))
# 
# res_summ_protein %>% 
#   ggplot(aes(x = max_BETA, y=max_PPI))+
#   geom_point()
# 
# res_summ_protein %>% 
#   filter(max_PPI > 0.9) %>% 
#   ggplot(aes(x = max_BETA))+
#   geom_histogram()

#-------------------------------------------------------------------------------
# simulate residuals of Y with the same correlation structure


# calculate residuals
Y_scaled = scale(Y_real, center = TRUE, scale = FALSE)
X_scaled = scale(X_real)

residuals = Y_scaled - X_scaled %*% obj_atlasqtl_realData$beta_vb 

# correlation matrix of residuals
residual_cor_mat = cor(residuals)
hist(as.vector(residual_cor_mat))

# correlation matrix of Y
Y_real_corr = cor(Y_real)
hist(as.vector(Y_real_corr))
saveRDS(Y_real_corr, "/rds/user/yl2021/hpc-work/hotspot_sim/Y_real_corr.rds")

#save data
saveRDS(residual_cor_mat, "/rds/user/yl2021/hpc-work/hotspot_sim/residual_cor_mat.rds")
saveRDS(residuals, "/rds/user/yl2021/hpc-work/hotspot_sim/residuals.rds")
################################################################################
#-------------------------------------------------------------------------------
# Data simulation
# setting up parameters, these are the real scopes, but it will take quite some time to run 
# you can play with small ones first as you like


active_ratio_p = 0.005
active_ratio_q = 0.2
max_tot_pve = 0.25
sh2 = 5

p = 1000
q = ncol(Y_real)
n = row(X_real)


p = 1000
q = 2919
n = nrow(X_real)

list = simulate_withRealGenome(as.matrix(X_real), 
                               q=2919, 
                               protein_ls = colnames(Y)[3:ncol(Y)],
                               active_protein = protein_rank[1:qt],
                               residual_cor_mat = Y_real_corr,
                               active_ratio_p = 0.005, 
                               max_tot_pve = 0.25,
                               sh2 = 5,
                               seed = 2020)

#save
saveRDS(list, file.path(save_path, "sim_list.rds"))

#load
list = readRDS(file.path(save_path, "sim_list_3.rds"))

X_sim = as.matrix(list$X)
Y_sim = as.matrix(list$Y)
pat_sim = list$pat
beta_sim = list$beta



#check distribution of simulated BETA
#abs
hist(abs(beta_sim[which(beta_sim!=0, arr.ind = T)]), breaks = 20)
#org
hist(beta_sim[which(beta_sim!=0, arr.ind = T)], breaks = 50)
#check distribution of simulated correlations of Y
Y_sim_corr = cor(Y_sim)
hist(as.vector(Y_sim_corr))



#-------------------------------------------------------------------------------
# Run atlasQTL on simulated data
# the prior expectation and variance of the number of predictors associated with each response.
obj_atlasqtl =  atlasqtl(Y = Y_sim, X = X_sim,
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=1000,
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

#save
saveRDS(obj_atlasqtl, "/rds/user/yl2021/hpc-work/hotspot_sim/obj_atlasqtl_sim.rds")
#load
obj_atlasqtl = readRDS("/rds/user/yl2021/hpc-work/hotspot_sim/obj_atlasqtl_sim_3.rds")

#plot results

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

