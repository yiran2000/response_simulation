# generate replicate data using echoseq and checkout whether the simulated SNP data makes sense 
library(dplyr)
library(data.table)
library(ggplot2)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(patchwork)
library(reshape2)
source("data_simulation_tools/simulate_withRealGenome.R")

#Read in X
geno_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD"
geno = fread(file.path(geno_path, "imputed_chr1_ld1_filtered.raw"))
X = geno[ ,7:ncol(geno)] #randomly subset n number of samples
colnames(X) = sapply(strsplit(colnames(X), "_"), `[`, 1)

# #Function to calculate MAF
# calculate_maf <- function(snp_column) {
#   # Count the frequency of each allele
#   allele_freq <- table(round(snp_column)) / length(snp_column)
#   
#   # Minor allele frequency is the minimum allele frequency
#   maf <- min(allele_freq)
#   return(maf)
# }
# 
# # Calculate MAF for each SNP
# maf_X <- apply(X, 2, mean)*0.5
# hist(maf_X)

#Replicate X by adding in noise and round to 0
add_noise = function(col, err=0.001){
  round(col + rnorm(length(col), mean=0, sd=err))
}
X_new = apply(X[1:35000, 1:1000], 2, add_noise)
colnames(X_new) = paste0("SNP_",1:ncol(X_new))
saveRDS(X_new, "replicated_SNPs_new.rds")

#Make LD plot of X_new
ld_matrix <- cor(X_new)
ld_melted <- reshape2::melt(ld_matrix)
ld_plot <- ggplot(ld_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(-1, 1), space = "Lab", name="LD (r)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "SNPs", y = "SNPs", title = "LD Plot")


#Replicate X
X_replicate = replicate_real_snps(n=35000, real_snps=as.matrix(X[1:35000, 1:1000]), bl_lgth=1000, p = NULL, maf_thres = NULL)
SNP_replicate = X_replicate[[1]]
maf_replicate = X_replicate[[2]]

#Make LD plot of X
ld_matrix <- cor(X[1:35000, 1:1000])
ld_melted <- reshape2::melt(ld_matrix)
ld_plot <- ggplot(ld_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(-1, 1), space = "Lab", name="LD (r)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "SNPs", y = "SNPs", title = "LD Plot")



#Make LD plot of replicated X
ld_matrix2 <- cor(SNP_replicate)
ld_melted2 <- reshape2::melt(ld_matrix2)
ld_plot2 <- ggplot(ld_melted2, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(-1, 1), space = "Lab", name="LD (r)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "SNPs", y = "SNPs", title = "LD Plot")



#Simulate based on the real SNPs and obtain running results
list = simulate_withRealGenome(X[1:35000, 1:1000], n = 35000, p=1000, q=3000, 
                               missing_ratio = 0,
                               active_ratio_q = 0.02, 
                               active_ratio_p = 0.05, 
                               max_tot_pve = 0.3,
                               sh2 = 20,
                               seed = NULL)
X_reformed = list$X
Y = list$Y
pat = list$pat

colnames(SNP_replicate) = colnames(X_reformed) 
res_atlas_X = atlasqtl(as.matrix(Y), as.matrix(X_reformed),
                       p0 = c(1, 4),
                       user_seed = 1, maxit= 50000,
                       thinned_elbo_eval = F,
                       batch = "y",
                       tol = 0.1,
                       # anneal_tol = c(1, 0.1),
                       anneal_tol = NULL,
                       # anneal = c(1, 2, 10),
                       anneal = NULL,
                       burn_in = 20,
                       epsilon_lb = c(2, 1.5, 0.25),
                       epsilon_it = c(0.1, 50, 0), #k, x_0, m
                       partial_elbo = F,
                       eval_perform = T)

#Simulate based on the replicated SNPs and obtain running results
list2 = simulate_withRealGenome(X_new, n = 35000, p=1000, q=3000, 
                               missing_ratio = 0,
                               active_ratio_q = 0.02, 
                               active_ratio_p = 0.05, 
                               max_tot_pve = 0.3,
                               sh2 = 20,
                               seed = NULL)
SNP_replicate_reformed = list2$X
Y = list2$Y
pat = list2$pat
res_atlas_X_replicate = atlasqtl(as.matrix(Y), as.matrix(SNP_replicate_reformed),
                       p0 = c(1, 4),
                       user_seed = 1, maxit= 50000,
                       thinned_elbo_eval = F,
                       batch = "y",
                       tol = 0.1,
                       # anneal_tol = c(1, 0.1),
                       anneal_tol = NULL,
                       # anneal = c(1, 2, 10),
                       anneal = NULL,
                       burn_in = 20,
                       epsilon_lb = c(2, 1.5, 0.25),
                       epsilon_it = c(0.1, 50, 0), #k, x_0, m
                       partial_elbo = F,
                       eval_perform = T)

saveRDS(res_atlas_X, "res_atlas_X.rds")
saveRDS(res_atlas_X, "res_atlas_X_replicate.rds")


p1 = melt(pat2) %>% group_by(Var1) %>% 
  summarise(n_assoc = sum(value == T)) %>% 
  ggplot(aes(x = Var1, y = n_assoc))+
  geom_point()

#Manhattan plot of atlasQTL results
gam_vb = res_atlas_X_replicate$gam_vb
rownames(gam_vb) = 1:nrow(gam_vb)
p2 = reshape2::melt(gam_vb) %>% 
  group_by(Var1) %>% 
  summarise(n_assoc = sum(value > 0.99)) %>% 
  ggplot(aes(x = Var1, y = n_assoc))+
  geom_point()
p1+p2
