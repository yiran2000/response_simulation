save_path = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/Sim3"
summary_ls = gtools::mixedsort(list.files(file.path(save_path, "summary")))
summary_df = lapply(summary_ls, function(f) read.table(file.path(save_path, "summary", f))) %>% 
  rbindlist(fill=TRUE)


summary_df = summary_df %>% 
  mutate(precision = TP/(FP+TP), 
         recall = TP/(TP+FN), 
         FPR = FP/(FP+TN)) 
  # left_join(param_table %>% mutate(task_id = 1:nrow(param_table)), by = c("data_id", "rep_id", "task_id","p","active_ratio_q")) %>% 
  # drop_na()


metrics <- c("iterations", "runtime", "ELBO", "AUROC","AUPRC","memory_peak",
             "TP","FP","TN","FN",
             "precision","recall", "FPR")


df_full <- summary_df %>%
  dplyr::filter(method == "Vanilla") %>% 
  dplyr::select(-task_id, -method, -iterations_full) %>%
  rename_with(~ paste0(.x, "_Vanilla"), all_of(metrics)) 

df_other <- summary_df%>% 
  dplyr::filter(method != "Vanilla") %>%
  dplyr::select(-task_id, -iterations_full) %>%
  rename_with(~ paste0(.x, "_other"), all_of(metrics)) 

# Join full metrics to the other methods by seed
df_combined <- df_other  %>% 
  left_join(df_full, by = c("data_id", "hm", "q", "active_ratio_q")) %>% drop_na()



# 
# # Separate full and non-full rows
# metrics <- c("iterations", "runtime", "ELBO", "AUROC","AUPRC","memory_peak")
# df_full <- summary_df %>%
#   filter(method == "full") %>% 
#   dplyr::select(-task_id, -method,-iterations_full) %>% 
#   rename_with(~ paste0(.x, "_full"), all_of(metrics)) 
# 
# df_other <- summary_df %>% 
#   # filter(method != "full") %>% 
#   filter(method == "partial_eIter") %>% 
#   dplyr::select(-task_id, -iterations_full) %>% 
#   rename_with(~ paste0(.x, "_other"), all_of(metrics)) 
# 
# # Join full metrics to the other methods by seed
# df_combined <- df_other  %>% 
#   left_join(df_full, by = c("rep_id","active_ratio_q","n","p","q")) %>% drop_na()



###############################################################################
#Speed comparison
p = df_combined %>%
  ggplot(aes(
    x = active_ratio_q,
    y = -(runtime_other/iterations_other - runtime_Vanilla/iterations_Vanilla) / (runtime_Vanilla/iterations_Vanilla),
    group = active_ratio_q
  )) +
  geom_boxplot() +
  scale_x_continuous(breaks = unique(summary_df$active_ratio_q))+
  ylab("Percentage runtime gain")+
  theme_bw()

ggsave(plot = p, filename = file.path(save_path, "plots/res_multi_q.pdf"),
       device = "pdf", 
       width = 5, height =4, dpi = 1200)

df_combined %>%
  ggplot(aes(
    x = active_ratio_q,
    y = (iterations_other - iterations_Vanilla) / iterations_Vanilla,
    group = active_ratio_q
  )) +
  geom_boxplot() +
  scale_x_continuous(breaks = unique(summary_df$active_ratio_q))+
  ylab("Percentage running time gain")+
  theme_bw()


#precision
p2 = summary_df %>% 
  ggplot(aes(x = active_ratio_q, y = precision, color = method))+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  theme_bw()

#recall
p3 = summary_df %>% 
  # filter(active_ratio_q == "0.005") %>% 
  ggplot(aes(x = active_ratio_q, y = recall, color = method))+
  facet_wrap(~p+hm, labeller = label_both)+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  theme_bw()




p1 + p2 + p3 +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom") 

#-------------------------------------------------------------------------------
#Speed gain w.r.t. active_ratio_q, per iterations




p2 = df_combined %>%
  ggplot(aes(
    x = active_ratio_q,
    y = (precision_other - precision_full) / (precision_full),
    color = as.factor(p)
  )) +
  ylab("Percentage difference in precision")+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  facet_wrap(~hm, labeller = label_both)+
  theme_bw()

p3 = df_combined %>%
  ggplot(aes(
    x = active_ratio_q,
    y = (recall_other - recall_full) / (recall_full),
    color = as.factor(p)
  )) +
  ylab("Percentage difference in recall")+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  facet_wrap(~hm, labeller = label_both)+
  theme_bw()


################################################################################
#runtime, relative difference:
p1 = df_combined %>%
  # filter(hm != 0.05) %>% 
  ggplot(aes(
    x = active_ratio_q,
    y = -(runtime_other - runtime_full) / (runtime_full)
  )) +
  ylab("Percentage running time gain")+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  facet_wrap(~hm, labeller = label_both)+
  theme_bw()

p2 = df_combined %>%
  # filter(hm != 0.05) %>% 
  ggplot(aes(
    x = active_ratio_q,
    y = -(iterations_other - iterations_full) / (iterations_full)
  )) +
  ylab("Percentage reduction of iterations")+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  facet_wrap(~hm, labeller = label_both)+
  theme_bw()

p3 = df_combined %>%
  # filter(hm != 0.05) %>% 
  ggplot(aes(
    x = active_ratio_q,
    y = -(runtime_other/iterations_other - runtime_full/iterations_full) / (runtime_full/iterations_full)
  )) +
  ylab("Percentage runtime again per iteration")+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  facet_wrap(~hm, labeller = label_both)+
  theme_bw()

#precision
p4 = summary_df %>% 
  # filter(hm != 0.05) %>% 
  ggplot(aes(x = active_ratio_q, y = precision, color = method))+
  facet_wrap(~hm, labeller = label_both)+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  theme_bw()

#recall
p5 = summary_df %>% 
  # filter(hm != 0.05) %>% 
  ggplot(aes(x = active_ratio_q, y = recall, color = method))+
  facet_wrap(~hm, labeller = label_both)+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  theme_bw()

#fpr
p6 = summary_df %>% 
  # filter(hm != 0.05) %>% 
  ggplot(aes(x = active_ratio_q, y = FPR, color = method))+
  facet_wrap(~hm, labeller = label_both)+
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar") +
  theme_bw()

p = p1 + p2 + p3 + p4+ p5+p6+
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "bottom") 

ggsave(plot = p, filename = file.path(save_path, "plots/res_multi.pdf"),
       device = "pdf", 
       width = 10, height = 5, dpi = 1200)
