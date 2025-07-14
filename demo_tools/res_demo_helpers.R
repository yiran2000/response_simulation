require(PRROC)
require(ggplot2)


#evaluate AUROC and AUPRC
ROC.atlasqtl = function(res_atlas, pat_sim){
  roc <- PRROC::roc.curve(scores.class0 = as.vector(res_atlas$gam_vb), weights.class0 = as.numeric(as.vector(as.matrix(pat_sim))), curve =TRUE)
  AUROC = roc$auc
  
  pr <- PRROC::pr.curve(scores.class0 = as.vector(res_atlas$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
  AUPRC = pr$auc.integral
  
  return(list(
    AUROC = AUROC,
    AUPRC = AUPRC
  ))
}

reshape_fittingCurve = function(fittingList, pat_sim = NULL){
  df = do.call(rbind, fittingList) %>% as.data.frame()
  df$iter <- 1:nrow(df) 
  df_long <- pivot_longer(df, -iter, names_to = "series", values_to = "value") %>% 
    mutate(associated = if_else(series %in% paste0("V",which(colSums(pat_sim) >0)), T, F))
  
  return(df_long)
}


resReshape.atlasqtl = function(obj_atlasqtl, X_sim, snp_gene_info){
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

#manhattan plots
manhattan.simData = function(pat_sim, snp_gene_info){
  p = data.frame(ID = rownames(pat_sim), ProteinCount = rowSums(pat_sim)) %>% 
    left_join(snp_gene_info, by = "ID") %>% 
    ggplot(
      aes(x = POS, y = ProteinCount)
    ) +
    geom_point()+
    theme_bw()
  
  return(p)
}

#enter the reshaped atlasqtl list
manhattan.atlasqtl = function(list_alas){
  p1 = list_alas[[1]]%>% 
    group_by(POS, ID) %>% 
    summarise(n_assoc = sum(corr_metric > 0.9)) %>% 
    ggplot(aes(x = POS, y = n_assoc))+
    geom_point()+
    theme_bw()
  
  p2 = list_alas[[2]] %>% 
    ggplot(aes(x = POS, y = theta))+
    geom_point()+
    theme_bw()
  
  return(p1, p2)
} 


#theta similarity, dot product
cos_sim <- function(u, v) {
  # Ensure dimensions match
  if (!all(dim(u) == dim(v))) {
    stop("u and v must have the same dimensions.")
  }

  dot_product <- u * v
  norm_u <- sqrt(u^2)
  norm_v <- sqrt(v^2)
  return(dot_product / (norm_u * norm_v))

}


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


