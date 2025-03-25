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


# simulation_tool_path = "~/project/response_simulation"
simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source(file.path(simulation_tool_path, "reshape_atlasQTL_res.R"))

demo_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/atlas_res_demo"
source(file.path(demo_tool_path, "demo_helpers/res_demo_helpers.R"))



################################################################################
#read in pre-save data
#X
list = readRDS(file.path(simulation_tool_path, "data/sim_list_3.rds"))
X_real = as.matrix(list$X)

#Y_real_corr
Y_real_corr = readRDS(file.path(simulation_tool_path, "data/Y_real_corr.rds"))

#protein rank
# protein_rank = readLines("data/protein_rank.txt")

#snp gene info
snp_gene_info = fread(file.path(simulation_tool_path, "data/snp_gene_info.csv"))

################################################################################
#simulate one single hotspot
list_sim = simulate_withRealGenome(
  X = X_real[,801:1000], 
  q = 1000,
  protein_ls = paste0("protein_", 1:1000),
  active_protein = NULL,
  active_ratio_p = 2/200, 
  active_ratio_q = 0.05,
  residual_cor_mat = NULL,
  vec_rho_phenos = NULL,
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  200, 
  min_cor = 0, max_cor = 0.9,
  m = 0.02,
  sh2 = 5,
  seed  = 2024
)


X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat

# Y_sim_corr = cor(Y_sim)
# hist(as.vector(beta_sim)[as.vector(beta_sim)!=0])

#visualize simulated associations
data.frame(ID = rownames(pat_sim), ProteinCount = rowSums(pat_sim)) %>% 
  left_join(snp_gene_info, by = "ID") %>% 
  ggplot(
    aes(x = POS, y = ProteinCount)
  ) +
  geom_point()+
  theme_bw()+
  ggtitle("Simulated # of associated proteins by SNP")


################################################################################
#run atlasQTL

# p = ncol(X_sim)
# q = 2919
# n = nrow(X_sim)


mu_t = 1
v_t = 4

devtools::load_all()

time_full = system.time({
  res_atlas_full = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                               p0 = c(mu_t, v_t),
                               user_seed = 1, maxit=1000,
                               batch = "y",
                               tol = 0.01, 
                               max_partial = 1,
                               end_full = T,
                               epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                               eval_perform = T,
                               thinned_elbo_eval = T)
})
#64.14s

system.time({
  res_atlas_partial = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                               p0 = c(mu_t, v_t),
                               user_seed = 1, maxit=1000,
                               batch = "y",
                               tol = 0.1, 
                               max_partial = Inf,
                               end_full = F,
                               epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                               eval_perform = T,
                               thinned_elbo_eval = T)
})



#plot(res_atlas)


################################################################################
# iterations
res_atlas_full$it
res_atlas_partial$it

# AUROC, AUPRC
ROC.atlasqtl(res_atlas_full, pat_sim)
ROC.atlasqtl(res_atlas_partial, pat_sim)


# fitting curve
res_atlas_full$perform_df %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial))+
  geom_point()+theme_bw()

res_atlas_partial$perform_df %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial))+
  geom_point()+theme_bw()

#final ELBO
res_atlas_full$lb_opt - res_atlas_partial$lb_opt

#theta similarity, dot product
cos_sim <- function(u, v) {
  # Ensure dimensions match
  if (!all(dim(u) == dim(v))) {
    stop("u and v must have the same dimensions.")
  }
  
  # If both are vectors
  if (is.vector(u) && is.vector(v)) {
    return(sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
  }
  
  # If both are matrices
  if (is.matrix(u) && is.matrix(v)) {
    # Row-wise cosine similarity
    dot_product <- rowSums(u * v)
    norm_u <- sqrt(rowSums(u^2))
    norm_v <- sqrt(rowSums(v^2))
    return(dot_product / (norm_u * norm_v))
  }
  
  stop("u and v must be both vectors or both matrices.")
}


cos_sim(res_atlas_full$theta_vb, res_atlas_partial$theta_vb) #because they locate on different local maxima, I suppose
cos_sim(res_atlas_full$zeta_vb, res_atlas_partial$zeta_vb)
cos_sim(res_atlas_full$gam_vb, res_atlas_partial$gam_vb) %>% hist()

################################################################################
#Theta curve full v.s. partial
x = seq(0,5, 0.001)
y = lognormal_cdf(x, mu=2, sigma=1.5, m=0.25)
plot(x,y)









################################################################################
#FDR
# exact FDR
gam_vb = obj_atlasqtl$gam_vb
predicted <- gam_vb > 0.9

TP <- sum((predicted==1) & (pat_sim == 1))
FP <- sum(predicted & (pat_sim == 0))
FDR <- FP / sum(predicted==1)
FDR

power = TP/sum(pat_sim == 1)
power


# loci-wise FDR
active_set_pred = which(apply(gam_vb, 2, function(col) any(col > 0.9)))
active_set_true = which(apply(pat_sim, 2, function(col) any(col == 1)))
TP_proteins = sum(active_set_pred %in% active_set_true)
FP_proteins = sum(!active_set_pred %in% active_set_true)

FP_proteins / (FP_proteins + TP_proteins)

TP_proteins/1000

