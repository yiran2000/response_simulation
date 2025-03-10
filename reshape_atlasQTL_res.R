reshape_atlasQTL_res = function(obj_atlasqtl, X_sim, snp_gene_info){
  # beta
  sigma_X = data.frame(
    ID = names(apply(X_sim, 2, sd)),
    sd = apply(X_sim, 2, sd))
  
  beta_mat = as.data.frame(obj_atlasqtl$beta_vb)
  beta_mat$ID = rownames(beta_mat)
  beta_df = reshape2::melt(beta_mat, id.vars = "ID", variable.name = "protein_name", value.name = "BETA")%>% 
    left_join(sigma_X, by = "ID") %>% 
    mutate(BETA_rescal = BETA/sd)
  
  
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
