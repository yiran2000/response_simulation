library(dplyr)
library(data.table)
library(ggplot2)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(patchwork)

#change the file directory as in your case
geno_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD"

source("data_simulation_tools/simulate_withRealGenome.R")
#-------------------------------------------------------------------------------
# setting up parameters, these are the real scopes, but it will take quite some time to run 
# you can play with small ones first as you like
q = 3000
p = 1000
n = 35000


active_ratio_p = 0.2
max_tot_pve = 0.3
sh2 = 20

#-------------------------------------------------------------------------------
#Data simulation
#X
geno = fread(file.path(geno_path, "imputed_chr1_ld1_filtered.raw"))
X = geno[ ,7:ncol(geno)] #randomly subset n number of samples
colnames(X) = sapply(strsplit(colnames(X), "_"), `[`, 1)

#with the seed as p, the dataset simulated should be the each for each p
list = simulate_withRealGenome(X, n = nrow(X), p=p, q=q, 
                               missing_ratio = 0,
                               active_ratio_q = 0.02, #sugges to fix this one
                               active_ratio_p = active_ratio_p, 
                               max_tot_pve = max_tot_pve,
                               sh2 = sh2,
                               seed = 2024)
X = list$X
Y = list$Y
pat = list$pat


# Mean impute Y
# If you want, you can set missing_ratio = 0.2 and then mean-impute, which is more simuler to the real case
# Y_meanImp = na_mean(as.data.frame(Y))


# manhattan plot
# Plot the number of proteins associated with each SNP to understand the simulated data better
# You can plot the one inferred by atlasQTL to compare
data.frame(SNP = 1:nrow(pat), ProteinCount = rowSums(pat)) %>% ggplot(
  aes(x = SNP, y = ProteinCount)
) +geom_point()



#-------------------------------------------------------------------------------
# Run atlasQTL
# the prior expectation and variance of the number of predictors associated with each response.
mu_t = 1
v_t = 4

# system.time(obj_atlasqtl <- atlasqtl(as.matrix(Y), as.matrix(X),
#                                   p0 = c(mu_t, v_t),
#                                   user_seed = 1, maxit= 50000,
#                                   thinned_elbo_eval = F,
#                                   batch = "y",
#                                   tol = 0.1,
#                                   # anneal_tol = c(1, 0.1),
#                                   anneal_tol = NULL,
#                                   # anneal = c(1, 2, 10),
#                                   anneal = NULL,
#                                   burn_in = 20,
#                                   epsilon_lb = c(2, 1.5, 0.25),
#                                   epsilon_it = c(0.1, 50, 0), #k, x_0, m
#                                   partial_elbo = F,
#                                   eval_perform = T))

obj_atlasqtl =  atlasqtl(as.matrix(Y_discovery[, 3:ncol(Y_discovery)]), as.matrix(X_discovery[, 2:ncol(X_discovery)]),
                                p0 = c(mu_t, v_t),
                                user_seed = 1, maxit=5000,
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



#-------------------------------------------------------------------------------
#Evaluate and save result
roc <- PRROC::roc.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat)), curve =TRUE)
AUROC = roc$auc

pr <- PRROC::pr.curve(scores.class0 = as.vector(obj_atlasqtl$gam_vb), weights.class0 = as.vector(as.matrix(pat)), curve =TRUE)
AUPRC = pr$auc.integral

output = 
  data.table(
    X_chunk_size = X_chunk_size,
    max_tot_pve = max_tot_pve,
    sh2 = sh2,
    # partial_update = partial_update,
    seed = seed,
    n_iter = obj_atlasqtl$it,
    ELBO = obj_atlasqtl$lb_opt,
    runtime_min = runtime[3]/60,
    AUROC = AUROC, 
    AUPRC = AUPRC,
    ave_assoc_per_snp = mean(rowSums(pat)),
    ave_assoc_per_protein = mean(colSums(pat)),
    ave_assoc_per_active_snp = mean(rowSums(pat)[rowSums(pat)>0]),
    ave_assoc_per_active_protein = mean(colSums(pat)[colSums(pat)>0])
  )

write.csv(output, file.path(save_path, paste0("output_", task_id, '.csv')), row.names = FALSE)
# saveRDS(obj_atlasqtl, file.path(save_path, paste0("obj_atlasqtl_", task_id, '.rds')))

