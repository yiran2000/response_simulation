#Simulate with only one active SNP and extremely low h^2, and large amount of associated proteins
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

simulation_tool_path = "~/project/response_simulation"
source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source(file.path(simulation_tool_path, "demo_tools/res_demo_helpers.R"))

list = readRDS(file.path("/rds/user/yl2021/hpc-work/hotspot_sim", "sim_list_3.rds"))
X_real = as.matrix(list$X)

p = 200
q = 1000

list_sim = simulate_withRealGenome(
  X = X_real[,1:p], 
  q = q,
  protein_ls = paste0("protein_", 1:q),
  active_protein = NULL,
  active_ratio_p = 1/p, 
  active_ratio_q = 200/q,
  residual_cor_mat = NULL,
  vec_rho_phenos = NULL,
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  50, 
  min_cor = 0, max_cor = 0.5,
  beta1 = 1,
  hm = 0.0005,
  seed  = 1
)

X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat


# saveRDS(list_sim, "/rds/user/yl2021/hpc-work/hotspot_sim/list_sim_example.rds")

res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                     p0 = c(1,4),
                     user_seed = 25, maxit=3000,
                     batch = "y",
                     tol = tol, 
                     start_partial = Inf, #Inf, 1000, 500, 100
                     max_partial = 1,
                     end_full = F, 
                     subsample_scheme = "geometric",
                     min_epsilon = 0,
                     geom_alpha = 0.97,
                     eval_perform = T,
                     thinned_elbo_eval = F)


#-------------------------------------------------------------------------------
res_atlas = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                     p0 = c(1,4),
                     user_seed = 25, maxit=3000,
                     batch = "y",
                     tol = tol, 
                     start_partial = Inf, #Inf, 1000, 500, 100
                     max_partial = 1,
                     end_full = F, 
                     subsample_scheme = "geometric",
                     min_epsilon = 0,
                     geom_alpha = 0.97,
                     eval_perform = T,
                     thinned_elbo_eval = F)


saveRDS(res_atlas, "/rds/user/yl2021/hpc-work/hotspot_sim/sim_res_atlas.rds")
atlas_df = process_corr_mat(obj_atlasqtl = res_atlas,
                            CHR = 19, protein_names =protein_names, snp_gene_info=snp_gene_info, protein_gene_info=protein_gene_info) 

                                                                                                                                                                                                          #On one simulated data, compare GWAS and atlasQTL
#GWAS
#-------------------------------------------------------------------------------
#Save simulated Y_sim into pheno files to run PLINK GWAS
Y = readRDS("/rds/user/yl2021/hpc-work/myukbb/proteomics_clean/residuals_meanImp_discovery.rds") %>% as.data.table()
ID = fread("/rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample")
pheno_sim = as.data.table(Y_sim) %>% 
  mutate(FID = Y$eid, IID = Y$eid) %>% 
  dplyr::select(FID, IID, everything())
pheno_sim_reordered = pheno_sim[match(ID$ID_1, pheno_sim$IID), ]

write.table(pheno_sim_reordered, file.path("/rds/user/yl2021/hpc-work/myukbb/my_pheno/pheno_sim.fam"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)

#save snplist
snp_ls = colnames(X_real[,1:p])
write.table(snp_ls, file = "/rds/user/yl2021/hpc-work/hotspot_sim/simulation_snplist.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#run GWAS
system("rm -rf /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim*")
system("find /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/ -type f -delete")
plink_cmd_GWAS <- paste(
  "plink2",
  "--bgen /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.bgen ref-first",
  "--sample /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample",
  "--extract  /rds/user/yl2021/hpc-work/hotspot_sim/simulation_snplist.txt",
  "--chr 19",
  "--pheno /rds/user/yl2021/hpc-work/myukbb/my_pheno/pheno_sim.fam",
  "--out /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim",
  "--glm"
)

system(plink_cmd_GWAS)


#load GWAS res
gwas_res = read_PLINK_GWAS("/rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19", pattern = "chr19_sim\\.PHENO\\d+\\.glm\\.linear", protein_names=colnames(Y_sim))
gwas_mt = as.matrix(gwas_res[[1]])
gwas_mt[gwas_mt == 0] = 1e-99

# beta
beta_mat = as.data.frame(gwas_res[[2]])
beta_mat$ID = rownames(beta_mat)
beta_df = reshape2::melt(beta_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")

# pval
corr_mat = as.data.frame(gwas_mt)
corr_mat$ID = rownames(corr_mat)

# joint
gwas_df = reshape2::melt(corr_mat, id.vars = "ID", variable.name = "protein_name", value.name = "corr_metric")%>%
  #Combine snp_gene_info
  left_join(snp_gene_info, by = "ID") %>%
  left_join(beta_df, by = c("ID", "protein_name")) %>%
  mutate(corr_metric = -log10(corr_metric))



#-------------------------------------------------------------------------------
# exact precision, recall, FPR with threshold
threshold = -log10(5e-8/1000)
labels <- as.numeric(as.vector(as.matrix(pat_sim)))  # ground truth (0 or 1)

predicted_col <- as.numeric(apply(-log10(gwas_mt) > threshold, 2, any))  #GWAS
predicted_col <- as.numeric(apply(res_atlas$gam_vb > 0.2, 2, any))  #alasQTL

labels_col     <- as.numeric(apply(pat_sim == 1, 2, any))

TP_col <- sum(predicted_col == labels_col & labels_col == 1)
TN_col <- sum(predicted_col == labels_col & labels_col == 0)
FN_col <- sum(predicted_col != labels_col & labels_col == 1)
FP_col <- sum(predicted_col != labels_col & labels_col == 0)

precision_col = TP_col /(FP_col +TP_col )
recall_col  = TP_col /(TP_col +FN_col)
FPR_col  = FP_col /(FP_col +TN_col)

#-------------------------------------------------------------------------------
#plot the auroc curve
labels_col     <- as.numeric(apply(pat_sim == 1, 2, any))
scores_col_gwas = (apply(-log10(gwas_mt), 2, max))
scores_col_atlas = apply(res_atlas$gam_vb, 2, max)
# Compute ROC curves
roc_gwas <- roc.curve(scores.class0 = scores_col_gwas[labels_col == 1],
                      scores.class1 = scores_col_gwas[labels_col == 0],
                      curve = TRUE)

roc_atlas <- roc.curve(scores.class0 = scores_col_atlas[labels_col == 1],
                       scores.class1 = scores_col_atlas[labels_col == 0],
                       curve = TRUE)

# Combine into a data.frame for ggplot
df_gwas <- data.frame(Specificity = roc_gwas$curve[,1],
                      Sensitivity = roc_gwas$curve[,2],
                      Method = "glm")

df_atlas <- data.frame(Specificity = roc_atlas$curve[,1],
                       Sensitivity = roc_atlas$curve[,2],
                       Method = "atlasQTL")

df_combined <- rbind(df_gwas, df_atlas)

# AUC values for annotation
auc_gwas <- round(roc_gwas$auc, 3)
auc_atlas <- round(roc_atlas$auc, 3)

# Plot
p1 = ggplot(df_combined %>% mutate(Method = recode(Method, glm = "PLINK")), aes(x = Specificity, y = Sensitivity, color = Method)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  annotate("text", x = 0.3, y = 0.15,
           label = paste0("AUC (PLINK): ", auc_gwas, "\nAUC (atlasQTL): ", auc_atlas),
           hjust = 0, size =3.5) +
  coord_fixed() +  # Ensure square aspect ratio
  scale_color_manual(values = c("PLINK" = "#0072B2", "atlasQTL" = "#D55E00")) +
  labs(title = "ROC Curve",
       x = "1 - Specificity", y = "Sensitivity",
       color = "Method") +
  theme_bw(base_size = 12) +
  ylim(0, 1) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 11), 
    legend.title = element_text(size = 12)
  )

#-------------------------------------------------------------------------------
#plot the auprc curve
# Compute PR curves
pr_gwas <- pr.curve(scores.class0 = scores_col_gwas[labels_col == 1],
                    scores.class1 = scores_col_gwas[labels_col == 0],
                    curve = TRUE)

pr_atlas <- pr.curve(scores.class0 = scores_col_atlas[labels_col == 1],
                     scores.class1 = scores_col_atlas[labels_col == 0],
                     curve = TRUE)

# Convert curves to data frames
df_gwas <- data.frame(Recall = pr_gwas$curve[, 1],
                      Precision = pr_gwas$curve[, 2],
                      Method = "glm")

df_atlas <- data.frame(Recall = pr_atlas$curve[, 1],
                       Precision = pr_atlas$curve[, 2],
                       Method = "atlasQTL")

df_combined <- rbind(df_gwas, df_atlas)

# AUC values for annotation
auc_gwas <- round(pr_gwas$auc.integral, 3)
auc_atlas <- round(pr_atlas$auc.integral, 3)

# Plot
p2 = ggplot(df_combined%>% mutate(Method = recode(Method, glm = "PLINK")), aes(x = Recall, y = Precision, color = Method)) +
  geom_line(size = 1) +
  annotate("text", x = 0.05, y = 0.15,
           label = paste0("AUC (PLINK): ", auc_gwas, "\nAUC (atlasQTL): ", auc_atlas),
           hjust = 0, size = 3.5) +
  coord_fixed() +  # Makes plot square
  scale_color_manual(values = c("PLINK" = "#0072B2", "atlasQTL" = "#D55E00")) +
  labs(title = "Precision-Recall Curve",
       x = "Recall", y = "Precision",
       color = "Method") +
  theme_bw(base_size = 12) +
  ylim(0, 1) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 11), 
    legend.title = element_text(size = 12)
  )


p <- ggarrange(
  p1, p2,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)


ggsave(plot = p, filename = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/resToy.jpg",
       device = "jpg", 
       width = 6.5, height = 4, dpi = 350)

# scores <- -log10(as.vector(gwas_mt))  # predicted scores
# predicted <- ifelse(scores > threshold, 1, 0)
# 
# #atlasQTL
# scores <- as.vector(res_atlas$gam_vb)  # predicted scores
# predicted <- ifelse(scores > 0.5, 1, 0)


# TP <- sum(predicted == labels & labels == 1)
# TN <- sum(predicted == labels & labels == 0)
# FN <- sum(predicted != labels & labels == 1)
# FP <- sum(predicted != labels & labels == 0)

# precision = TP/(FP+TP)
# recall = TP/(TP+FN)
# FPR = FP/(FP+TN)


beta_sim_mat= as.data.frame(beta_sim)
beta_sim_mat$ID = rownames(beta_sim_mat)
beta_sim_df = reshape2::melt(beta_sim_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")
p0 <- beta_sim_df %>%
  filter(BETA != 0) %>%
  ggplot(aes(x = BETA)) +
  geom_histogram(fill = "steelblue", color = "white", alpha = 0.5) +
  ggtitle(expression("Distribution of simulated " ~ beta)) +
  xlab(expression("Regression coefficient" ~ beta ~ "(non-zero)")) +
  theme_bw()

ggsave(plot = p0, filename = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/betaToy.pdf",
       device = "pdf", 
       width = 4, height = 4, dpi = 1200)

#-------------------------------------------------------------------------------
#ppi v.s. simulated beta
atlas_df %>% left_join(
  as.data.frame(beta_sim) %>% 
    mutate(ID = rownames(pat_sim)) %>% 
    reshape2::melt(id.vars = "ID", variable.name = "protein_name", value.name = "BETA_sim")%>%
    left_join(snp_gene_info, by = "ID"), 
  by = c("ID", "protein_name")
) %>% 
  ggplot(aes(x = corr_metric, y = BETA_sim)) +
  geom_point()


#atlas inferred beta v.s. simulated beta
atlas_df %>% left_join(
  gwas_df,
  by = c("ID", "protein_name")
) %>% 
  ggplot(aes(x = BETA.x, y = BETA.y, color = corr_metric.x)) +
  geom_point()


#PPI v.s. p-value
joint_df = atlas_df %>% left_join(
  gwas_df, 
  by = c("ID", "protein_name")
)%>% group_by(protein_name) %>% 
  summarise(
    corr_metric.x.max = max(corr_metric.x),
    corr_metric.y.max = max(corr_metric.y),
    BETA.x.max = max(abs(BETA.x)),
    BETA.y.max = max(abs(BETA.y))
  ) %>% left_join(
   data.frame(protein_name = colnames(pat_sim), 
              sig_sim = labels_col),
    by = c("protein_name")
  ) 


joint_df%>% 
  ggplot(aes(x = corr_metric.x.max, y = corr_metric.y.max, color = as.factor(sig_sim))) +
  geom_point(alpha = 0.7)+theme_bw()


joint_df%>% 
  ggplot(aes(x = BETA.x.max, y = BETA.y.max, color = as.factor(sig_sim))) +
  geom_point()

#-------------------------------------------------------------------------------
#manhattan plot of simulated v.s. GWAS inferred


p1 = as.data.frame(pat_sim) %>% 
  mutate(ID = rownames(pat_sim)) %>% 
  reshape2::melt(id.vars = "ID", variable.name = "protein_name", value.name = "if_assoc")%>%
  left_join(snp_gene_info, by = "ID") %>%
  dplyr::group_by(ID, POS) %>% 
  dplyr::summarise(
    n_assoc = sum(if_assoc)
  ) %>% 
  ggplot(aes(x = POS))+
  geom_segment(aes(xend = POS, yend = 0, y = n_assoc), color = "black", linewidth = 0.3, alpha = 0.5)+
  theme_bw()+
  ggtitle("Simulated associations")+
  xlab(NULL)+
  ylab("Nb. associated proteins")

p2 = gwas_df  %>% 
  dplyr::group_by(ID, POS) %>% 
  dplyr::summarise(
    n_assoc = sum(corr_metric >= -log10(5e-8/1000))
  ) %>% 
  ggplot(aes(x = POS))+
  geom_segment(aes(xend = POS, yend = 0, y = n_assoc), color = "black", linewidth = 0.3, alpha = 0.5)+
  theme_bw()+
  ggtitle("Univariate screening inferred associations")+
  xlab(NULL)+
  ylab("Nb. associated proteins")

#summarise by loci
loci = (gwas_df  %>% 
  dplyr::group_by(ID, POS) %>% 
  dplyr::summarise(
    n_assoc = sum(corr_metric >= -log10(5e-8/1000))
  ) %>% arrange(desc(n_assoc)) %>% pull(POS))[1]

gwas_df %>%
  filter(corr_metric >= threshold) %>%
  distinct(protein_name) %>% nrow()



p3 = atlas_df %>% 
  dplyr::group_by(ID, gene_POS) %>% 
  dplyr::summarise(
    n_assoc = sum(corr_metric >= 0.5)
  ) %>% 
  ggplot(aes(x = gene_POS))+
  geom_segment(aes(xend = gene_POS, yend = 0, y = n_assoc), color = "black", linewidth = 0.3, alpha = 0.5)+
  theme_bw()+
  ggtitle("atlasQTL inferred associations")+
  xlab(NULL)+
  ylab("Nb. associated proteins")

p1 + p2 + p3 + plot_layout(ncol = 1)



#p-value v.s. simulated beta
gwas_df %>% left_join(
  as.data.frame(beta_sim) %>% 
    mutate(ID = rownames(pat_sim)) %>% 
    reshape2::melt(id.vars = "ID", variable.name = "protein_name", value.name = "BETA_sim")%>%
    left_join(snp_gene_info, by = "ID"), 
  by = c("ID", "protein_name", "POS")
) %>% 
  ggplot(aes(x = corr_metric, y = BETA_sim)) +
  geom_point()

#GWAS inferred beta v.s. simulated beta
gwas_df %>% left_join(
  as.data.frame(beta_sim) %>% 
    mutate(ID = rownames(pat_sim)) %>% 
    reshape2::melt(id.vars = "ID", variable.name = "protein_name", value.name = "BETA_sim")%>%
    left_join(snp_gene_info, by = "ID"), 
  by = c("ID", "protein_name", "POS")
) %>% 
  ggplot(aes(x = BETA, y = BETA_sim, color = corr_metric)) +
  geom_point()





#-------------------------------------------------------------------------------
#LD plot

#first convert bgen to befile
system("plink2 --bgen /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.bgen ref-first --sample /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample --make-bed --out /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp")

# run the following from terminal:
# module unload plink
# module load plink/1.9
# plink \
# --bfile /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp \
# --extract /rds/user/yl2021/hpc-work/hotspot_sim/simulation_snplist.txt \
# --chr 19 \
# --r2 --ld-window 1000 --ld-window-kb 1000 --ld-window-r2 0 \
# --out /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp_ld


library(snpStats)
library(LDheatmap)

ld_results <- read.table("/rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp_ld.ld", header = TRUE)

ld_results = ld_results %>% filter(SNP_A %in% snp_ls, SNP_B %in% snp_ls)

# function that convert ld df to matrix
ld_convert = function(ld_results){
  snp_positions <- unique(c(ld_results$BP_A, ld_results$BP_B))  # Unique positions
  snp_positions <- sort(snp_positions)  # Sort positions
  
  # Initialize an empty matrix
  ld_matrix <- matrix(NA, nrow = length(snp_positions), ncol = length(snp_positions),
                      dimnames = list(snp_positions, snp_positions))
  
  # Fill the LD matrix with R2 values
  for (i in 1:nrow(ld_results)) {
    pos_a <- ld_results$BP_A[i]
    pos_b <- ld_results$BP_B[i]
    ld_matrix[as.character(pos_a), as.character(pos_b)] <- ld_results$R2[i]
    ld_matrix[as.character(pos_b), as.character(pos_a)] <- ld_results$R2[i]  # Symmetric
  }
  
  # Replace NA with 0 for visualization
  diag(ld_matrix) = 1
  
  return(ld_matrix)
}

ld_matrix = ld_convert(ld_results)

rgb.palette <- colorRampPalette(rev(c("#F0F0F0", "black")), space = "rgb")
LDheatmap(ld_matrix,
              genetic.distances = as.numeric(rownames(ld_matrix)),
              color=rgb.palette(18))






