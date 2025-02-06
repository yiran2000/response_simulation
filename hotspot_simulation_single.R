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
source("data_simulation_tools/simulate_withRealGenome.R")
source("reshape_atlasQTL_res.R")

################################################################################
#read in pre-save data
#X
list = readRDS("data/sim_list_3.rds")
X_real = as.matrix(list$X)

#Y_real_corr
Y_real_corr = readRDS("data/Y_real_corr.rds")

#protein rank
protein_rank = readLines("data/protein_rank.txt")

#snp gene info
snp_gene_info = fread("data/snp_gene_info.csv")

################################################################################
#simulate one single hotspot
list_sim = simulate_withRealGenome(X = X_real[,501:1000], q = 2919,
                                   protein_ls = colnames(Y_real_corr),
                                   active_SNP = NULL,
                                   active_protein = protein_rank[1:1000],
                                   active_ratio_p = 1/500, 
                                   active_ratio_q = NULL,
                                   residual_cor_mat = Y_real_corr,
                                   missing_ratio = 0, 
                                   max_tot_pve = NULL, 
                                   m = 0.001,
                                   # second shape parameter for the beta distribution controlling the hotspot propensities
                                   sh2 = 5,
                                   seed  = 2024)


X_sim = list_sim$X
Y_sim = list_sim$Y
beta_sim = list_sim$beta
pat_sim = list_sim$pat

# Y_sim_corr = cor(Y_sim)
hist(as.vector(beta_sim)[as.vector(beta_sim)!=0])

################################################################################
#run atlasQTL

# p = ncol(X_sim)
# q = 2919
# n = nrow(X_sim)


mu_t = mean(colSums(pat_sim))
v_t = mu_t*5

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

################################################################################

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
res_simData = reshape_atlasQTL_res(obj_atlasqtl,snp_gene_info = snp_gene_info)
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

