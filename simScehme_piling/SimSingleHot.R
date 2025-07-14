


################################################################################
#Using existing correlation in 2919 proteins
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


################################################################################
#equvicorrelated, self-controlled correlation
# active_ratio_per_block = 0.2
# active_protein = paste0("protein_", c(1:50, 101:150, 201:250, 301:350, 401:450))


list_sim = simulate_withRealGenome(
  X = X_real, 
  q = 1000,
  protein_ls = paste0("protein_", 1:1000),
  active_protein = paste0("protein_", c(1:100, 201:300, 401:500, 601:700, 801:900)),
  # active_protein = paste0("protein_", c(1, 201, 401, 601, 801)),
  active_ratio_p = 1/200, 
  active_ratio_q = NULL,
  residual_cor_mat = NULL,
  vec_rho_phenos = c(0.05, 0.25, 0.45, 0.65, 0.85),
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  200, 
  # min_cor = 0, max_cor = 0.9,
  hm = 0.001,
  sh2 = 5,
  seed  = 2024
)


# saveRDS(list_sim, "/rds/user/yl2021/hpc-work/hotspot_sim/list_sim_5.rds")
# list_sim = readRDS("/rds/user/yl2021/hpc-work/hotspot_sim/list_sim_4.rds")

X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat

beta_sim_mat= as.data.frame(beta_sim)
beta_sim_mat$ID = rownames(beta_sim_mat)
beta_sim_df = reshape2::melt(beta_sim_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")


active_protein = names(which(colSums(pat_sim)>0))
################################################################################
#Visualize correlation, simulated beta, etc.
#-------------------------------------------------------------------------------
#BETA
hist(as.vector(beta_sim)[as.vector(beta_sim)!=0])
max_beta = max(as.vector(beta_sim))

beta_sim_df %>% 
  filter(BETA !=0) %>% 
  ggplot(aes(x = BETA))+
  geom_histogram(fill = "steelblue", color = "white", alpha = 0.5) +
  ggtitle("Distribution of simulated BETA (non-zero)")+
  theme_bw()

#-------------------------------------------------------------------------------
#correlation in simulated Y
Y_sim_corr = cor(Y_sim)
data.frame(corr = as.vector(Y_sim_corr[upper.tri(Y_sim_corr)])) %>%
  filter(corr > 0) %>% 
  ggplot(aes(x = corr))+
  geom_histogram(fill = "steelblue", color = "white", alpha = 0.5) +
  ggtitle("Distribution of correlations in simulated Y")+
  theme_bw()


melted_Y_sim_corr <- reshape2::melt(as.matrix(Y_sim_corr))

ggplot(melted_Y_sim_corr, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Correlation Structure in simulated responses", fill = "Correlation")




#get the maximum correlation df
corr_df = as.data.frame(as.table(as.matrix(Y_sim_corr))) %>%
  rename(Variable1 = Var1, Variable2 = Var2, Correlation = Freq) %>%
  filter(Variable1 != Variable2)

active_protein = names(which(colSums(pat_sim) > 0))
inactive_protein = setdiff(paste0("protein_", 1:1000), active_protein)

block_assignments <- tibble(
  protein = colnames(pat_sim),
  block = ceiling(as.numeric(gsub("protein_", "", colnames(pat_sim))) / 200)
)

filtered_corr_df <- corr_df %>%
  left_join(block_assignments, by = c("Variable1" = "protein")) %>%
  rename(Block1 = block) %>%
  left_join(block_assignments, by = c("Variable2" = "protein")) %>%
  rename(Block2 = block) %>%
  filter(Block1 == Block2) %>%  # Ensure same block
  filter((Variable1 %in% active_protein & Variable2 %in% inactive_protein) |
           (Variable2 %in% active_protein & Variable1 %in% inactive_protein))


max_corr_df = filtered_corr_df %>% 
  group_by(Variable1) %>% 
  summarise(max_corr = max(Correlation)) %>% 
  rename(protein_name = "Variable1")
