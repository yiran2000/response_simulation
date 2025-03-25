#Visualize LD and correlation structure of proteins 


#LD
library(LDheatmap)
LD_matrix <- cor(X_sim, use = "pairwise.complete.obs")^2
LDheatmap(LD_matrix, genetic.distances = NULL, title = "LD Heatmap of SNPs")  


#protein
dist_matrix <- as.dist(1 - Y_real_corr)  # Convert correlation to a distance
hc <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering


ordered_matrix <- Y_real_corr[hc$order, hc$order]

mat = ordered_matrix[1:1000, 1:1000]

melted_data <- reshape2::melt(ordered_matrix, na.rm = TRUE)


# Plot using ggplot2

melted_data %>% 
  filter(Var1 %in% protein_rank[1:1000], Var2 %in% protein_rank[1:1000]) %>% 
  ggplot(aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Upper Triangular Heatmap of Correlation Blocks",
       x = "Variables", y = "Variables", fill = "Correlation")
