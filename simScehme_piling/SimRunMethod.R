################################################################################
#atlasQTL
mu_t =1
v_t = 4

# devtools::load_all()
obj_atlasqtl =  atlasqtl(Y = Y_sim, X = X_sim,
                         p0 = c(mu_t, v_t),
                         user_seed = 1, 
                         maxit = 3000,
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
gam_vb = obj_atlasqtl$gam_vb
res_simData = reshape_atlasQTL_res(obj_atlasqtl, X_sim,snp_gene_info)


################################################################################
#GWAS
#-------------------------------------------------------------------------------
#Save simulated Y_sim into pheno files to run PLINK GWAS
Y = readRDS("/rds/user/yl2021/hpc-work/myukbb/proteomics_clean/residuals_meanImp_discovery.rds") %>% as.data.table()
ID = fread("/rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample")
pheno_sim = as.data.table(Y_sim) %>% 
  mutate(FID = Y$eid, IID = Y$eid) %>% 
  dplyr::select(FID, IID, everything())
pheno_sim_reordered = pheno_sim[match(ID$ID_1, pheno_sim$IID), ]

write.table(pheno_sim_reordered, file.path(mybiobank, "my_pheno/pheno_sim.fam"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)

#run GWAS
system("rm -rf /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim*")
system("find /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/ -type f -delete")
plink_cmd_GWAS <- paste(
  "plink2",
  "--bgen /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.bgen ref-first",
  "--sample /rds/user/yl2021/hpc-work/hotspot_sim/simulated_snp.sample",
  "--chr 19",
  "--pheno /rds/user/yl2021/hpc-work/myukbb/my_pheno/pheno_sim.fam",
  "--out /rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19/chr19_sim",
  "--glm"
)

system(plink_cmd_GWAS)

#-------------------------------------------------------------------------------
#load GWAS res
gwas_res = read_PLINK_GWAS("/rds/user/yl2021/hpc-work/myukbb/UKB_GWAS/chr19", pattern = "chr19_sim\\.PHENO\\d+\\.glm\\.linear", protein_names=colnames(Y_sim))

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
