run_path = "/rds/user/yl2021/hpc-work/UKB_run_results/chr19_selected_anneal"
save_path = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme"
chunk_info$chunk_id = 1:nrow(chunk_info)
#Read in output files
output = lapply( mixedsort(list.files(run_path,pattern = "\\.csv$")), function(file){
  fread(file.path(run_path, file))
}) %>% rbindlist()

output


################################################################################
#Hotspot sizes
size_df = sun_pQTL %>% 
  group_by(`Region ID`) %>% 
  summarise(n = n()) %>% arrange(desc(n))

quantile(size_df$n)
hist(size_df$n)
  
################################################################################
#demonstrate similarlity in performance

ukb_res_files = gtools::mixedsort(list.files(run_path,pattern = "^obj_atlasqtl_.*\\.rds$"))
gam_vb_comb_full = lapply(1:6, function(i) {
  gam_vb = process_corr_mat(atlasqtl_path = file.path(run_path, ukb_res_files[i]),
                            CHR = CHR, protein_names =protein_names, snp_gene_info=snp_gene_info, protein_gene_info=protein_gene_info) %>% filter(corr_metric > 0.5)
  print(paste("chunk:", i, "processed"))
  return(gam_vb)
}) %>% rbindlist()


gam_vb_comb_partial = lapply(7:12, function(i) {
  gam_vb = process_corr_mat(atlasqtl_path = file.path(run_path, ukb_res_files[i]),
                            CHR = CHR, protein_names =protein_names, snp_gene_info=snp_gene_info, protein_gene_info=protein_gene_info) %>% filter(corr_metric > 0.5)
  print(paste("chunk:", i, "processed"))
  return(gam_vb)
}) %>% rbindlist()

# gam_vb_comb %>%
#   sample_frac(0.1) %>%
#   ggplot(aes(x = corr_metric.full, y = corr_metric.partial))+
#   geom_point()
# 
# sum (gam_vb_comb$corr_metric.full * gam_vb_comb$corr_metric.partial) / sum(sqrt(gam_vb_comb$corr_metric.full^2)*sqrt(gam_vb_comb$corr_metric.partial^2))

################################################################################
#loci-wise summarization

summarise_locus <- function(gam_vb, dist = 0.5e6, corr_thres = 0.5) {
  # Filter and rank by significance
  gam_vb$assigned <- FALSE
  gam_vb_filt <- gam_vb %>% filter(corr_metric > corr_thres) %>% arrange(desc(corr_metric))
  
  loci_list <- list()
  locus_id <- 1
  
  while (any(!gam_vb_filt$assigned)) {
    
    top_row <- gam_vb_filt[which(!gam_vb_filt$assigned)[1], ]
    center_pos <- top_row$gene_POS
    
    # Define locus as ±0.5Mb around center SNP
    in_locus <- with(gam_vb_filt, abs(gene_POS - center_pos) <= dist & !assigned)
    
    locus_range <- range(gam_vb_filt$gene_POS[in_locus])
    
    # Mark SNPs as assigned
    gam_vb_filt$assigned[in_locus] <- TRUE
    
    # Get all proteins in the full table within locus window
    in_window <- gam_vb$gene_POS >= locus_range[1] & gam_vb$gene_POS <= locus_range[2]

    # For each protein, get max beta and corr_metric

      
    loci_summary_long = gam_vb[in_window, ] %>%
      group_by(protein_name) %>%
      arrange(desc(corr_metric), desc(abs(BETA))) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        locus_start = locus_range[1],
        locus_end = locus_range[2],
        locus_POS = (locus_range[1] + locus_range[2]) / 2,
        protein_POS = (start_position + end_position) / 2,
        cis_trans = if_else(abs(locus_POS - protein_POS) <= 1e6, "cis", "trans"),
        snp_ID = ID
      ) %>%
      select(
        protein_name, locus_start, locus_end, snp_ID,
        snp_POS = gene_POS, gene_POS = locus_POS, BETA, corr_metric,
        protein_start = start_position, protein_end = end_position, protein_POS, cis_trans, chromosome_name
      )


    loci_list[[locus_id]] <- loci_summary_long
    locus_id <- locus_id + 1
  }
  

  loci_summary_long <- bind_rows(
    Map(function(df, i) {
      df$ID <- paste0("Locus_", i)
      df
    }, loci_list, seq_along(loci_list))
  )
  return(loci_summary_long)
}



loci_partial = summarise_locus(gam_vb_comb_partial) 
loci_full = summarise_locus(gam_vb_comb_full)

#-------------------------------------------------------------------------------
#Locus summarising
#summarising table of the partial
loci_partial %>% group_by(ID,locus_start,locus_end,cis_trans) %>% 
  summarise(n_protein = dplyr::n())%>%
  dplyr::mutate(ID_num = as.integer(gsub("Locus_", "", ID))) %>%
  arrange(ID_num) %>%
  pivot_wider(
    id_cols = c(ID_num, locus_start, locus_end),
    names_from = cis_trans,
    values_from = n_protein,
    values_fill = 0,   # fill missing cis/trans with 0
    names_prefix = "n_"
  )%>% 
  kableExtra::kbl(format = "latex", booktabs = TRUE)

mytheme = theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),  # removes both x and y major grid lines
    panel.grid.minor = element_blank(),  # removes both x and y minor grid lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "none"
  )

p1 = loci_full %>%
  left_join(loci_partial, by = c("locus_start","locus_end", "protein_name","ID")) %>%
  ggplot(aes(x = BETA.x, y = BETA.y)) +
  geom_point() +
  xlab("BETA (Vanilla CAVI)") +
  ylab("BETA (AF-CAVI)") +
  coord_fixed() +
  mytheme

p2 = loci_full %>%
  left_join(loci_partial, by = c("locus_start","locus_end", "protein_name","ID")) %>%
  ggplot(aes(x = corr_metric.x, y = corr_metric.y)) +
  geom_point() +
  xlab("PPI (Vanilla CAVI)") +
  ylab("PPI (AFIO-CAVI)") +
  coord_fixed() +
  mytheme



# 移除子图中的 xlab
p1_noxlab <- p1 + xlab(NULL)
p2_noxlab <- p2 + xlab(NULL)

# 对齐 panel + axis tick
aligned <- align_plots(p1_noxlab, p2_noxlab, align = "hv", axis = "tblr")

# 组合面板
combined_plot <- plot_grid(plotlist = aligned, nrow = 1)

# 添加两个 xlabel，分别控制位置
p = ggdraw() +
  draw_plot(plot_grid(plotlist = aligned, nrow = 1), x = 0, y = 0.10, width = 1, height = 0.88) +
  draw_label("BETA (Vanilla CAVI)", x = 0.3, y = 0.1, hjust = 0.5, size = 12, vjust = 1) +
  draw_label("PPI (Vanilla CAVI)",  x = 0.8, y = 0.1, hjust = 0.5, size = 12, vjust = 1)



ggsave(plot = p, filename = file.path("/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/res_Real1.jpg"),
       device = "jpg", 
       width = 6, height = 3, dpi = 350)


  


#-------------------------------------------------------------------------------
#manhattan plot
loci_partial %>% 
  group_by(ID,gene_POS) %>% 
  summarise(n_assoc = n()) %>% 
  ggplot(aes(x = gene_POS))+
  geom_segment(aes(xend = gene_POS, yend = 0, y = n_assoc), color = "black", linewidth = 0.3, alpha = 0.5)+
  scale_color_distiller(direction = 1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "grey80", linewidth = 0.1))+
  my_theme+
  ggtitle("Number of proteins associated")+
  xlab(NULL)+
  ylab("count")

  
plot_correlation(corr_df = loci_full , corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 200, plotting_thres = 0.5,
                        legend_theme = NULL,
                        title = expression("PPI of SNP-protein pairs"))[[1]]
  



p1 = plot_correlation(corr_df = loci_full , corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 50, plotting_thres = 0.5,
                 legend_theme = NULL,
                 title = expression("PPI of SNP-protein pairs"))[[3]]
p2 = plot_correlation(corr_df = loci_partial, corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 50, plotting_thres = 0.9,
                      legend_theme = NULL,
                      title = expression("PPI of SNP-protein pairs"))[[3]]


p = ggarrange(
  p1+theme(
    legend.position = "bottom", 
    legend.box = "horizontal",  # Ensures items are in a row
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  ),
  p2+theme(
    legend.position = "bottom", 
    legend.box = "horizontal",  # Ensures items are in a row
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  ), 
  ncol = 2, align = "h",labels = c("a", "b"), widths = c(1, 1)
)

ggsave(plot = p, filename = file.path(save_path, "res_RealData2_bothAnneal.pdf"),
       device = "pdf", 
       width = 14, height = 9, dpi = 1200)


#summarising table

sum(loci_full$protein_name != loci_partial$protein_name)

loci_full %>% group_by(ID,locus_start,locus_end) %>% summarise(n_protein = dplyr::n()) %>% 
  left_join(loci_partial %>% group_by(ID,locus_start,locus_end) %>% summarise(n_protein = dplyr::n()),
            by = c("ID","locus_start","locus_end")) %>% 
  mutate(ID_num = as.integer(gsub("Locus_", "", ID))) %>%
  arrange(ID_num) %>%
  select(-ID_num) %>% 
  mutate(num_overlap = n_protein.x) %>% 
  kableExtra::kbl(format = "latex", booktabs = TRUE, 
      caption = "Summary of proteins per locus in Vanilla and AF-CAVI")




################################################################################
p1 = plot_correlation(corr_df = gam_vb_comb_full, corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 50, plotting_thres = 0.9,
                      legend_theme = NULL,
                      title = expression("PPI of SNP-protein pairs"))[[3]]
p2 = plot_correlation(corr_df = gam_vb_comb_partial, corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 50, plotting_thres = 0.9,
                      legend_theme = NULL,
                      title = expression("PPI of SNP-protein pairs"))[[3]]


p = ggarrange(
  p1+theme(
    legend.position = "bottom", 
    legend.box = "horizontal",  # Ensures items are in a row
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  ),
  p2+theme(
    legend.position = "bottom", 
    legend.box = "horizontal",  # Ensures items are in a row
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  ), 
  ncol = 2, align = "h",labels = c("a", "b"), widths = c(1, 1)
)

ggsave(plot = p, filename = file.path(save_path, "res_RealData2_bothAnneal.pdf"),
       device = "pdf", 
       width = 12, height = 7, dpi = 1200)


################################################################################
#Check a particular chunk

chunk_ls = c(3,8,17,28,30,37)

i = 2
start = chunk_info[[chunk_ls[i], "start"]]
end = chunk_info[[chunk_ls[i], "stop"]]

p1 = plot_correlation(corr_df = gam_vb_comb_full %>% filter(gene_POS > start & gene_POS < end), corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 50, plotting_thres = 0.9,
                      legend_theme = NULL,
                      title = expression("PPI of SNP-protein pairs"))[[3]]
p2 = plot_correlation(corr_df = gam_vb_comb_partial%>% filter(gene_POS > start & gene_POS < end), corr_thres = 0.5, point_alpha = 0.7, point_size = 1.5, count_cap = 50, plotting_thres = 0.9,
                      legend_theme = NULL,
                      title = expression("PPI of SNP-protein pairs"))[[3]]


ggarrange(
  p1+theme(
    legend.position = "bottom", 
    legend.box = "horizontal",  # Ensures items are in a row
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  ),
  p2+theme(
    legend.position = "bottom", 
    legend.box = "horizontal",  # Ensures items are in a row
    legend.spacing.x = unit(0.5, "cm")  # Adjust spacing between legend items
  ), 
  ncol = 2, align = "h",labels = c("a", "b"), widths = c(1, 1)
)
