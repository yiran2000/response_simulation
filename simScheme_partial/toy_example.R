library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(tictoc)
library(PRROC)

library(atlasqtl)
library(echoseq)

################################################################################
#load functions

HPC = F

if(HPC){
  simulation_tool_path = "~/project/response_simulation"
}else{
  simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
}

source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source(file.path(simulation_tool_path, "demo_tools/res_demo_helpers.R"))


#Load in Real SNP data
list = readRDS(file.path(simulation_tool_path, "data/sim_list_3.rds"))
X_real = as.matrix(list$X)


################################################################################
#Toy example: 200 SNPs, 1000 responses, 35000 samples, only 1 hotspot associated with 10 proteins

list_sim = simulate_withRealGenome(
  X = X_real[1:35000,801:1000], 
  q = 1000,
  protein_ls = paste0("protein_", 1:1000),
  active_protein = NULL,
  active_ratio_p = 1/200, 
  active_ratio_q = 10/1000,
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

saveRDS(list_sim, file.path(simulation_tool_path, "data/list_sim_toy.rds"))

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

time_partial = system.time({
  res_atlas_partial = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                               p0 = c(mu_t, v_t),
                               user_seed = 1, maxit=1000,
                               batch = "y",
                               tol = 0.01, 
                               max_partial = Inf,
                               end_full = F,
                               epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                               eval_perform = T,
                               thinned_elbo_eval = T)
})

time_partial_full = system.time({
  res_atlas_partial_full = atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                               p0 = c(mu_t, v_t),
                               user_seed = 1, maxit=1000,
                               batch = "y",
                               tol = 0.01, 
                               max_partial = Inf,
                               end_full = T,
                               epsilon = c(0.2, 1, 5, 1),#e0, emax, ec50, n
                               eval_perform = T,
                               thinned_elbo_eval = T)
})

saveRDS(res_atlas_full, file.path(simulation_tool_path, "data/res_toy_full.rds"))
saveRDS(res_atlas_partial, file.path(simulation_tool_path, "data/res_toy_partial.rds"))
saveRDS(res_atlas_partial_full, file.path(simulation_tool_path, "data/res_toy_partial_full.rds"))


#plot(res_atlas)

################################################################################
#Latex summarising table

summary_df <- data.frame(
  Method = c("Full", "Partial", "Partial fullEval"),
  Iterations = c(res_atlas_full$it, res_atlas_partial$it, res_atlas_partial_full$it),
  runtime = c(time_full[[3]], time_partial[[3]], time_partial_full[[3]]),
  ELBO = c(res_atlas_full$lb_opt, res_atlas_partial$lb_opt, res_atlas_partial_full$lb_opt),
  AUROC = c(ROC.atlasqtl(res_atlas_full, pat_sim)$AUROC, ROC.atlasqtl(res_atlas_partial, pat_sim)$AUROC, ROC.atlasqtl(res_atlas_partial_full, pat_sim)$AUROC),
  AUPRC = c(ROC.atlasqtl(res_atlas_full, pat_sim)$AUPRC, ROC.atlasqtl(res_atlas_partial, pat_sim)$AUPRC, ROC.atlasqtl(res_atlas_partial_full, pat_sim)$AUPRC),
  theta_similarty = c(NA, mean(cos_sim(res_atlas_full$theta_vb, res_atlas_partial$theta_vb)), mean(cos_sim(res_atlas_full$theta_vb, res_atlas_partial_full$theta_vb))),
  zeta_similarty = c(NA, mean(cos_sim(res_atlas_full$zeta_vb, res_atlas_partial$zeta_vb)), mean(cos_sim(res_atlas_full$zeta_vb, res_atlas_partial_full$zeta_vb))),
  gamma_similarity = c(NA, mean(cos_sim(res_atlas_full$gam_vb, res_atlas_partial$gam_vb)), mean(cos_sim(res_atlas_full$gam_vb, res_atlas_partial_full$gam_vb)))
)


library(knitr)
kable(summary_df, format = "latex", caption = "Summary of Full, Partial, Partial fullEval version of the CAVI algorithm", booktabs = TRUE)



#-------------------------------------------------------------------------------
# visualize cosine similarity of gamma

# # Calculate similarities
# sim1 <- cos_sim(res_atlas_full$gam_vb, res_atlas_partial$gam_vb)
# sim2 <- cos_sim(res_atlas_full$gam_vb, res_atlas_partial_full$gam_vb)
# 
# # Find common range
# x_min <- min(min(sim1), min(sim2))
# x_max <- max(max(sim1), max(sim2))
# bin_width <- (x_max - x_min) / 30
# 
# # Plot 1
# p1 <- data.frame(similarity = sim1) %>% 
#   ggplot(aes(x = similarity)) +
#   geom_histogram(fill = "#619CFF", color = "white", binwidth = bin_width) +
#   theme_bw() +
#   labs(title = "Distribution of Cosine Similarity (Partial)",
#        x = "Cosine Similarity",
#        y = "Frequency") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_x_continuous(limits = c(x_min, x_max))
# 
# # Plot 2
# p2 <- data.frame(similarity = sim2) %>% 
#   ggplot(aes(x = similarity)) +
#   geom_histogram(fill = "#619CFF", color = "white", binwidth = bin_width) +
#   theme_bw() +
#   labs(title = "Distribution of Cosine Similarity (Partial Full)",
#        x = "Cosine Similarity",
#        y = "Frequency") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_x_continuous(limits = c(x_min, x_max))


################################################################################


################################################################################
#Runnig course visualizations:

# ELBO
p1.1 = res_atlas_full$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  # mutate(partial = factor("full", levels = c("full", "partial"))) %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial))+
  labs(x = "Iteration",
       y = "ELBO",
       title = "Full")+
  geom_point()+theme_bw()

p1.2 = res_atlas_partial$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial))+
  labs(x = "Iteration",
       y = "ELBO",
       title = "Partial")+
  geom_point()+theme_bw()

p1.3 = res_atlas_partial_full$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial))+
  labs(x = "Iteration",
       y = "ELBO",
       title = "Partial fullEval")+
  geom_point()+theme_bw()


p1 <-  wrap_plots(p1.1, p1.2, p1.3, guides = "collect") +   theme(legend.position = "bottom")
  # plot_annotation(tag_levels = 'a') &
  

# ELBO diff
p2.1 = res_atlas_full$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = ELBO_diff, color = partial))+
  labs(x = "Iteration",
       y = "ELBO difference",
       title = "Full")+
  geom_point()+
  theme_bw()

p2.2 = res_atlas_partial$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = ELBO_diff, color = partial))+
  labs(x = "Iteration",
       y = "ELBO difference",
       title = "Partial")+
  geom_point()+
  theme_bw()

p2.3 = res_atlas_partial_full$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = ELBO_diff, color = partial))+
  labs(x = "Iteration",
       y = "ELBO difference",
       title = "Partial fullEval")+
  geom_point()+
  theme_bw()

p2 <- wrap_plots(p2.1, p2.2, p2.3, guides = "collect") +theme(legend.position = "bottom")

p = ggarrange(
  p1,
  p2, 
  ncol = 1, align = "v",labels = c("a", "b")
)
ggsave(file.path(simulation_tool_path, "simScheme_partial/plots/ELBO_compare.png"),p, width = 12, height = 6)






#Time per loop and Eplison of the partial algorithm
p3.1 = res_atlas_partial$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = e))+
  labs(x = "Iteration",
       y = "Epsilon")+
  geom_point()+
  theme_bw()
p3.2 = res_atlas_partial$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = time_loop))+
  labs(x = "Iteration",
       y = "Runtime per loop (s)")+
  geom_point()+
  theme_bw()
p3.3 = res_atlas_partial$perform_df %>% 
  mutate(partial = if_else(subsample_size == 1000, F, T)) %>% 
  ggplot(aes(x = iter, y = subsample_size))+
  labs(x = "Iteration",
       y = "Subsample size")+
  geom_point()+
  theme_bw()

p3 = ggarrange(
  p3.1,
  p3.2, 
  p3.3,
  nrow = 1, align = "hv",labels = c("a", "b","c"),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(file.path(simulation_tool_path, "simScheme_partial/plots/partial_process.png"),p3, width = 9, height = 3)





