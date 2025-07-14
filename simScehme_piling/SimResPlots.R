################################################################################
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
  summarise(n_assoc = sum(corr_metric > -log10(0.05/(1000*1000)))) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  geom_point()+
  theme_bw()+
  ggtitle("Inferred # of associated proteins by SNP (GWAS)")

#theta
p1.3 = res_simData[[2]] %>%
  ggplot(aes(x = POS, y = theta))+
  geom_point()+
  theme_bw()+
  ggtitle("Inferred theta")



p1 = wrap_plots(p1.1, p1.2, p1.3, ncol = 1)


################################################################################
#Checking piling up effect of atlasQTL

pat_df = data.frame(
  protein_name = names(colSums(pat_sim)>0),
  if_active = colSums(pat_sim)>0
)

atlas_summ = 
  res_simData[[1]] %>% 
  left_join(beta_sim_df, by = c("ID", "protein_name")) %>%
  group_by(protein_name) %>% 
  drop_na() %>% 
  summarise(max_BETA.x = max(abs(BETA_rescal)),
            max_BETA.y = max(abs(BETA.y)),
            max_ppi = max(corr_metric)) %>% 
  mutate(block = ceiling(as.numeric(gsub("protein_", "", colnames(pat_sim))) / 200)) %>% 
  mutate(corr = case_when(
    block == 1 ~ 0.05,
    block == 2 ~ 0.25,
    block == 3 ~ 0.45,
    block == 4 ~ 0.65,
    block == 5 ~ 0.85
  ))  %>% 
  left_join(pat_df, by = "protein_name") 
  
res_simData[[1]] %>% 
  left_join(beta_sim_df, by = c("ID", "protein_name")) %>% 
  ggplot(aes(x = BETA.x, y = BETA.y))+
  geom_point()


#simulated BETA v.s. inferred BETA
p1 = atlas_summ %>%
  ggplot(aes(x = abs(max_BETA.x), y = abs(max_BETA.y), color = max_ppi))+
  geom_point(alpha = 0.5)+
  facet_wrap(~corr, nrow = 1)+
  xlab("BETA (atlasQTL)")+
  ylab("BETA (simulated)")+
  theme_bw()


#simulated BETA v.s. inferred PPI
atlas_summ%>%
  ggplot(aes(x = max_ppi, y = max_BETA.y, color = if_active))+
  geom_point(alpha = 0.5)+
  xlab("PPI (atlasQTL)")+
  ylab("BETA (sim)")+
  theme_bw()

atlas_summ %>% 
  filter(protein_name %in% paste0("protein_", 1:100)) %>%
  ggplot(aes(x = max_ppi, y = abs(max_BETA.x), color = if_active))+
  geom_point(alpha = 0.5)+
  xlab("PPI (atlasQTL)")+
  ylab("BETA (simulated)")+
  theme_bw()

#Those that are highly correlated with the true proteins
#maximum correlation
# Assume Y_sim_corr is already your correlation matrix


p3 = atlas_summ %>% 
  ggplot(aes(x = max_ppi, y = corr, color = if_active))+
  geom_point(alpha = 0.5)+
  xlab("PPI (atlasQTL)")+
  ylab("correlation with active proteins")+
  theme_bw()

################################################################################
#GWAS 
gwas_summ = gwas_df %>% 
  left_join(beta_sim_df, by = c("ID", "protein_name")) %>%
  group_by(protein_name) %>% 
  drop_na() %>% 
  summarise(max_BETA.x = max(BETA.x),
            max_BETA.y = max(BETA.y),
            max_pval = max(corr_metric)) %>%
  mutate(block = ceiling(as.numeric(gsub("protein_", "", colnames(pat_sim))) / 200)) %>% 
  mutate(corr = case_when(
    block == 1 ~ 0.05,
    block == 2 ~ 0.25,
    block == 3 ~ 0.45,
    block == 4 ~ 0.65,
    block == 5 ~ 0.85
  ))  %>% 
  left_join(pat_df, by = "protein_name")


#GWAS beta sim v.s. inf
p2 = gwas_summ %>% 
  ggplot(aes(x = abs(max_BETA.x), y = abs(max_BETA.y), color = max_pval))+
  geom_point(alpha = 0.5)+
  facet_wrap(~corr, nrow = 1)+
  xlab("BETA (GWAS)")+
  ylab("BETA (simulated)")+
  theme_bw()


p4 = gwas_summ %>% 
  ggplot(aes(x = max_pval, y = corr, color = if_active))+
  geom_point(alpha = 0.5)+
  xlab("-log10(pval) (GWAS)")+
  ylab("correlation with active proteins")+
  theme_bw()


################################################################################
#atlasQTL v.s. GWAS
atlas_gwas_joint_df = res_simData[[1]] %>% 
  left_join(gwas_df, by = c("protein_name", "ID","POS"))

atlas_gwas_joint_summ =  atlas_gwas_joint_df %>% 
  group_by(protein_name) %>% 
  drop_na() %>% 
  summarise(max_ppi = max(corr_metric.x),
            max_logpval = max(corr_metric.y),
            max_BETA.x = max(abs(BETA_rescal)),
            max_BETA.y = max(abs(BETA.y))) %>% 
  left_join(pat_df, by = "protein_name") %>% 
  mutate(block = ceiling(as.numeric(gsub("protein_", "", colnames(pat_sim))) / 200)) %>% 
  mutate(corr = case_when(
    block == 1 ~ 0.05,
    block == 2 ~ 0.25,
    block == 3 ~ 0.45,
    block == 4 ~ 0.65,
    block == 5 ~ 0.85
  )) 


#PPI v.s. BETA
atlas_gwas_joint_summ %>%   
  ggplot(aes(x =max_ppi, y = max_logpval, color = if_active))+
  facet_wrap(~corr, nrow = 1)+
  xlab("maximum PPI by protein")+
  ylab(expression("maximum -log"[10]*"(p-value) by protein"))+
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = -log10(1.7e-8), linetype = "dashed", linewidth = 0.2)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2)+
  ggtitle("Maximum -log10pval v.s. PPI")+
  labs(color = "If simulated active")+
  theme_bw()


#BETA atlas and GWAS
atlas_gwas_joint_summ %>%   
  ggplot(aes(x =max_BETA.x, y = max_BETA.y, color = if_active))+
  facet_wrap(~corr, nrow = 1)+
  xlab("maximum PPI by protein")+
  ylab(expression("maximum -log"[10]*"(p-value) by protein"))+
  geom_point(alpha = 0.5)+
  ggtitle("Maximum -log10pval v.s. PPI")+
  labs(color = "If simulated active")+
  theme_bw()
