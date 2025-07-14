library(paletteer)
library(scales)
library(kableExtra)
library(ggrastr)
library(patchwork)


# param_table = param_table %>% rename(method = "scheme")

save_path = "/rds/user/yl2021/hpc-work/partial_update_validate_simScheme/Sim1_new_new"
summary_ls = gtools::mixedsort(list.files(file.path(save_path, "summary")))
summary_df = lapply(summary_ls, function(f) read.table(file.path(save_path, "summary", f))) %>% 
  rbindlist(fill=TRUE)


summary_df$method <- factor(summary_df$method, levels = c(
  "Vanilla",
  "RF",
  "AFE",
  "AFI",
  "AFIO"
))

summary_df = summary_df %>% 
  mutate(precision = TP/(FP+TP), 
         recall = TP/(TP+FN), 
         FPR = FP/(FP+TN))

#filter out outliers, keep those with all 5 methods finished
summary_df2 = summary_df %>% 
  group_by(task_id) %>%
  filter(all(c("Vanilla", "RF", "AFE", "AFI", "AFIO") %in% method),#filter out those with all 5 methods 
         !task_id %in% c(302,484,499)) 

# remove duplicates, and select 50 unique ones
task_selected <- summary_df2 %>%
  group_by(hm, active_ratio_q, data_id) %>%
  slice_max(task_id, n = 1, with_ties = FALSE) %>% #always take the larger task_id
  group_by(hm, active_ratio_q) %>%
  slice_sample(n = 50) %>%
  ungroup() %>% pull(task_id) %>% sort()

summary_df3 = summary_df2 %>%
  filter(task_id %in% task_selected)

summary_df =summary_df3 %>% ungroup()

#check
check = summary_df %>% group_by(hm, active_ratio_q,method) %>% summarise(n = n_distinct(data_id))

#check whether each of the hm, active_ratio_q group for different methods contain the same data_ids
# Step 1: Count data_id sets per group
data_id_check <- summary_df %>%
  group_by(hm, active_ratio_q, method) %>%
  summarise(
    data_ids = list(sort(unique(data_id))),
    n_data_id = n_distinct(data_id),
    .groups = "drop"
  )

# Step 2: Compare sets within each (hm, active_ratio_q)
data_id_check_consistency <- data_id_check %>%
  group_by(hm, active_ratio_q) %>%
  group_modify(~ {
    ref <- .x$data_ids[[1]]
    .x %>%
      mutate(is_equal = map_lgl(data_ids, ~ identical(.x, ref)))
  }) %>%
  ungroup()

data_id_check_consistency$is_equal


#--------------------------------------------------------------------------------
# comparision table
metrics <- c("iterations", "runtime", "ELBO", "AUROC","AUPRC","memory_peak",
             "TP","FP","TN","FN",
             "precision","recall", "FPR",
             "time_init", "time_loop", "time_local", "time_global")


df_full <- summary_df %>%
  filter(method == "Vanilla") %>%
  dplyr::select(-method) %>%
  rename_with(~ paste0(.x, "_Vanilla"), all_of(metrics))

df_other <- summary_df%>% filter(method != "Vanilla") %>%
  # dplyr::select(-task_id) %>%
  rename_with(~ paste0(.x, "_other"), all_of(metrics))


# Join full metrics to the other methods by seed
df_combined <- df_other  %>% 
  left_join(df_full, by = c("task_id", "data_id","hm", "active_ratio_q","tol"))

# #remove outliers
# task_outlier = df_combined %>% 
#   mutate(rel_loop = -(time_loop_other/iterations_other - time_loop_Vanilla/iterations_Vanilla) / (time_loop_Vanilla/iterations_Vanilla)) %>% 
#   filter(rel_loop < 0) %>% pull(task_id) %>% unique()

#check
check = df_combined %>% group_by(hm, active_ratio_q,method) %>% summarise(n = n_distinct(data_id))

#save
write.csv(summary_df, file.path(save_path, "summary_df.csv"), row.names = FALSE)
write.csv(df_combined, file.path(save_path, "df_combined.csv"), row.names = FALSE)


################################################################################
paletteer_d("colorBlindness::Blue2DarkOrange12Steps", 6)
custom_colors <- c(
  "Vanilla" = "black",    
  "RF" = "darkgrey",
  "AFE" = "#51C3CCFF",  
  "AFI" = "#FFAD65FF",         
  "AFIO" = "#FF8E32FF"
)

mytheme = theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),  # removes both x and y major grid lines
    panel.grid.minor = element_blank(),  # removes both x and y minor grid lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "none"
  )


################################################################################
#Time division in Vanilla CAVI
summary_df %>% 
  filter(method == "Vanilla") %>% 
 mutate(
  local_ratio = time_local/(runtime - time_init),
  global_ratio = time_global/(runtime - time_init),
  elbo_ratio = 1 - (time_global+time_local)/(runtime - time_init)
) %>%
  # group_by(active_ratio_q, hm) %>% 
  summarise(
    Local  = sprintf("%.4f $\\pm$ %.4f", mean(local_ratio, na.rm = TRUE), sd(local_ratio, na.rm = TRUE)),
    Global = sprintf("%.4f $\\pm$ %.4f", mean(global_ratio, na.rm = TRUE), sd(global_ratio, na.rm = TRUE)),
    ELBO   = sprintf("%.4f $\\pm$ %.4f", mean(elbo_ratio, na.rm = TRUE), sd(elbo_ratio, na.rm = TRUE))
  ) 

# %>%
#   kable(format = "latex", booktabs = TRUE, escape = FALSE,
#         caption = "Time ratio summary by active ratio $q$ (mean $\\pm$ std deviation).",
#         col.names = c("Local", "Global", "ELBO"))


################################################################################
#Main figure 1:
#hm = 0.15, relative time gain in local updates


p1 = df_combined %>%
  mutate(rel_loop = -(time_loop_other - time_loop_Vanilla) / (time_loop_Vanilla),
         a_q_label = factor(paste0("a[q]==", active_ratio_q))) %>% 
  filter(hm == 0.15) %>% 
  ggplot(aes(
    y = -100 * (time_local_other - time_local_Vanilla) / (time_local_Vanilla),
    x = method,
    color = method
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~a_q_label, labeller = label_parsed, nrow = 1)+
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors)+
  ylab(expression(Delta ~ "runtime (local) %"))+
  scale_y_continuous(breaks = seq(-50,100, 25))+
  mytheme

p2 = df_combined %>%
  mutate(rel_loop = -(time_loop_other - time_loop_Vanilla) / time_loop_Vanilla,
         a_q_label = factor(paste0("a[q]==", active_ratio_q))) %>% 
  filter(hm == 0.15) %>% 
  ggplot(aes(
    y = -100* ((runtime_other) - (runtime_Vanilla)) / (runtime_Vanilla),
    x = method,
    color = method
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~a_q_label, labeller = label_parsed, nrow = 1)+
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors)+
  ylab(expression(Delta ~ "runtime (total) %"))+
  scale_y_continuous(breaks = seq(-50,100, 25))+
  mytheme

p3 = df_combined %>%
  mutate(rel_loop = -(time_loop_other - time_loop_Vanilla) / time_loop_Vanilla,
         a_q_label = factor(paste0("a[q]==", active_ratio_q))) %>% 
  filter(hm == 0.15) %>% 
  ggplot(aes(
    y = 100*(iterations_other - iterations_Vanilla) / iterations_Vanilla,
    x = method,
    color = method
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~a_q_label, labeller = label_parsed, nrow = 1)+
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors)+
  ylab(expression(Delta ~ "Iterations %"))+
  scale_y_continuous(breaks = seq(-50,100, 25))+
  mytheme


p = (p1+theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank())+ 
     p2+theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
     p3)+
     plot_layout(nrow = 3)


ggsave(plot = p, filename = file.path(save_path, "plots/res_runtime.jpg"),
       device = "jpg", 
       width = 7, height = 6, dpi = 350)


#-------------------------------------------------------------------------------
#place full results in the appendix

metric_sel = c("time_local", "runtime", "iterations")

summary_long <- summary_df %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>% filter(metric %in% metric_sel)


#relative value
vanilla_values <- summary_long %>%
  filter(method == "Vanilla") %>%
  select(hm, active_ratio_q, data_id, metric, vanilla_value = value)

#calculate relative difference
summary_relative <- summary_long %>%
  left_join(vanilla_values, by = c("hm", "active_ratio_q", "data_id", "metric")) %>%
  mutate(rel_diff = round((value - vanilla_value) / vanilla_value,2))%>% drop_na()

#summary and format
summary_formatted = summary_relative %>%
  filter(metric %in% metric_sel, method != "Vanilla") %>%
  group_by(method,hm,active_ratio_q,metric) %>%
  summarise(metric_mean = mean(rel_diff*100),
            metric_sd = sd(rel_diff*100)) %>%
  mutate(metric_mean = if_else(metric == "runtime"| metric == "time_local", -metric_mean, metric_mean)) %>% 
  mutate(
    rel_summary = sprintf("%.2f ± %.2f", metric_mean,  metric_sd)
  ) 

#long to wide
method_order <- c("RF", "AFE", "AFI", "AFIO")

summ_runtime_full = summary_formatted %>%
  select(method, hm, active_ratio_q, metric, rel_summary) %>%
  pivot_wider(names_from = metric, values_from = rel_summary) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(hm, active_ratio_q, method) %>%  
  relocate(method, .after = active_ratio_q) 

# summary_print %>% 
#   kbl(format = "latex", booktabs = TRUE)
# 


################################################################################
#Table 1: ELBO, precision, recall, FPR relative difference and absolute value

metric_sel = c("ELBO", "precision","recall", "FPR")

summary_long <- summary_df %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>% filter(metric %in% metric_sel)

#absolute value, main context demo only
summary_long %>% 
  filter(hm == 0.15, active_ratio_q == 0.005, metric != "ELBO") %>% 
  group_by(method, metric) %>% 
  summarise(metric_mean = mean(value),
            metric_sd = sd(value))  %>%
  mutate(
    rel_summary = if_else(
      metric == "FPR",
      sprintf("%.2e ± %.2e", metric_mean, metric_sd),  
      sprintf("%.2f ± %.2f", metric_mean, metric_sd)
    )
  )%>% 
  select(method, metric, rel_summary) %>%
  pivot_wider(names_from = metric, values_from = rel_summary)%>% 
  kbl(format = "latex", booktabs = TRUE)


#for the appendix
summ_metric_full = summary_long %>% 
  filter(metric != "ELBO") %>% 
  group_by(active_ratio_q, hm, method, metric) %>% 
  summarise(metric_mean = mean(value),
            metric_sd = sd(value))  %>%
  mutate(
    rel_summary = if_else(
      metric == "FPR",
      sprintf("%.2e ± %.2e", metric_mean, metric_sd),  
      sprintf("%.2f ± %.2f", metric_mean, metric_sd)
    )
  )%>% 
  select(method, metric, rel_summary) %>%
  pivot_wider(names_from = metric, values_from = rel_summary)

summ_metric_full %>%
  left_join(summ_runtime_full, by = c("hm", "active_ratio_q", "method")) %>% 
  select(hm, active_ratio_q,method, iterations,runtime,time_local,FPR,precision,recall) %>% 
  arrange(hm, active_ratio_q)%>% 
  kbl(format = "latex", booktabs = TRUE) %>% 
  collapse_rows(columns = 1:2, 
                row_group_label_position = "identity",
                latex_hline = "major")  # adds lines between (active_ratio_q, hm)

#-------------------------------------------------------------------------------

# #relative value
# vanilla_values <- summary_long %>%
#   filter(method == "Vanilla") %>%
#   select(hm, active_ratio_q, data_id, metric, vanilla_value = value)
# 
# #calculate relative difference
# summary_relative <- summary_long %>%
#   left_join(vanilla_values, by = c("hm", "active_ratio_q", "data_id", "metric")) %>%
#   mutate(rel_diff = round((value - vanilla_value) / vanilla_value,2))%>% drop_na()
# 
# #summary and format
# summary_formatted = summary_relative %>%
#   filter(metric %in% metric_sel, method != "Vanilla") %>%
#   group_by(method,hm,active_ratio_q,metric) %>%
#   summarise(metric_mean = mean(rel_diff * 100),
#             metric_sd = sd(rel_diff * 100)) %>%
#   mutate(
#     rel_summary = sprintf("%.2f ± %.2f", metric_mean,  metric_sd)
#   ) 
# 
# #long to wide
# method_order <- c("RF", "AFE", "AFI", "AFIO")
# 
# summary_print = summary_formatted %>%
#   select(method, hm, active_ratio_q, metric, rel_summary) %>%
#   pivot_wider(names_from = metric, values_from = rel_summary) %>%
#   mutate(method = factor(method, levels = method_order)) %>%
#   arrange(hm, active_ratio_q, method) %>%  
#   relocate(method, .after = active_ratio_q) 


#-------------------------------------------------------------------------------
#plot 
p1 = df_combined %>%
  mutate(a_q_label = factor(paste0("a[q]==", active_ratio_q))) %>% 
  filter(hm == 0.15) %>% 
  ggplot(aes(
    y =  (precision_other - precision_Vanilla) / precision_Vanilla,
    x = method,
    color = method
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~a_q_label, labeller = label_parsed)+
  # scale_color_manual(values = custom_colors)+
  # scale_fill_manual(values = custom_colors)+
  ylab(expression(Delta ~ "Precision"))+
  mytheme

p2 = df_combined %>%
  mutate(rel_loop = -(time_loop_other - time_loop_Vanilla) / time_loop_Vanilla,
         a_q_label = factor(paste0("a[q]==", active_ratio_q))) %>% 
  filter(hm == 0.15) %>% 
  ggplot(aes(
    y =  (FPR_other - FPR_Vanilla) / FPR_Vanilla,
    x = method,
    color = method
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~a_q_label, labeller = label_parsed)+
  # scale_color_manual(values = custom_colors)+
  # scale_fill_manual(values = custom_colors)+
  ylab(expression(Delta ~ "FPR"))+
  mytheme

p3 = df_combined %>%
  mutate(rel_loop = -(time_loop_other - time_loop_Vanilla) / time_loop_Vanilla,
         a_q_label = factor(paste0("a[q]==", active_ratio_q))) %>% 
  filter(hm == 0.15) %>% 
  ggplot(aes(
    y =  (recall_other - recall_Vanilla) / recall_Vanilla,
    x = method,
    color = method
  )) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~a_q_label, labeller = label_parsed)+
  # scale_color_manual(values = custom_colors)+
  # scale_fill_manual(values = custom_colors)+
  ylab(expression(Delta ~ "Recall"))+
  mytheme

p1/p2/p3


################################################################################
#In one simulation:
param_table %>% 
  filter(active_ratio_q == 0.005, hm == 0.15)

perform_df_full = read.table(file.path(save_path,"perform/perform_186_Vanilla.csv"))
perform_df_AFIO = read.table(file.path(save_path,"perform/perform_186_AFIO.csv"))

p1 = perform_df_full %>% 
  filter(annealing == F) %>% 
  ggplot(aes(x = iter, y = time_local))+
  geom_point()+
  scale_y_continuous(limits = c(0,3.5))+
  theme_bw()

p2 = perform_df_AFIO %>% 
  filter(annealing == F) %>% 
  ggplot(aes(x = iter, y = time_local))+
  geom_point()+
  scale_y_continuous(limits = c(0,3.5))+
  theme_bw()

p3 = perform_df_AFIO %>% 
  filter(annealing == F) %>% 
  ggplot(aes(x = iter, y = subsample_size))+
  geom_point()+
  theme_bw()

p1+p2+p3

sum(perform_df_AFIO$time_loop)
sum(perform_df_full$time_loop)

perform_df_AFIO %>% 
  filter(annealing == F) %>% 
  ggplot(aes(x = iter, y = time_local))+
  geom_point()+
  scale_y_continuous(limits = c(0,3))+
  theme_bw()


################################################################################
#Distribution of zeta
zeta_ls= lapply(gtools::mixedsort(list.files(file.path(save_path, "zeta"))), function(f) read.table(file.path(save_path, "zeta", f)))

zeta_df <- Reduce(function(x, y) left_join(x, y, by = "V1"), zeta_ls[13:32])
zeta_df_long <- zeta_df %>%
  pivot_longer(
    cols = -V1,  # all columns except protein names
    names_to = "replicate",
    values_to = "value"
  ) %>% 
  dplyr::select(-replicate)


ggplot(zeta_df_long, aes(x = V1, y = value)) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "steelblue",
    size = 2
  ) +
  stat_summary(
    fun.data = mean_cl_normal,  # automatically computes 95% CI
    geom = "errorbar",
    width = 0.2,
    color = "steelblue"
  ) +
  theme_minimal() +
  labs(x = "Protein", y = "Zeta (mean ± 95% CI)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

################################################################################
#training trajectories
# zeta_nm_ls = gtools::mixedsort(list.files(file.path(save_path, "zeta_fittingCurve")))
# zeta_df_ls= lapply(gtools::mixedsort(list.files(file.path(save_path, "zeta_fittingCurve"))), function(f) read.table(file.path(save_path, "zeta_fittingCurve", f)))
# 

zeta_df_full = read.table(file.path(save_path,"zeta_fittingCurve/zeta_fittingCurve_1.csv"))
zeta_df_partial_eIter_adaptive_anneal = read.table(file.path(save_path,"zeta_fittingCurve/zeta_fittingCurve_201.csv"))

# Define color mapping for logical values
zeta_colors <- c("TRUE" = "black", "FALSE" = "grey")

# p1 plot
p1 <- zeta_df_full %>%
  # filter((associated == FALSE & iter %% 10 == 0) | associated == TRUE) %>% 
  ggplot(aes(x = iter, y = value, group = series)) +
  rasterize(geom_point(data = subset(zeta_df_full, associated == FALSE),
                       aes(color = associated), alpha = 0.3, size = 0.3), dpi = 300) +
  rasterize( geom_point(data = subset(zeta_df_full, associated == TRUE),
                        aes(color = associated), alpha = 0.3, size = 0.3), dpi = 300) +  
  scale_color_manual(values = zeta_colors) +
  ggtitle("Vanilla CAVI") +
  theme_bw()

p2 <- zeta_df_partial_eIter_adaptive_anneal %>%
  ggplot(aes(x = iter, y = value, group = series)) +
  rasterize(geom_point(data = subset(zeta_df_full, associated == FALSE),
                       aes(color = associated), alpha = 0.3, size = 0.3), dpi = 300) +
  rasterize( geom_point(data = subset(zeta_df_full, associated == TRUE),
                        aes(color = associated), alpha = 0.3, size = 0.3), dpi = 300) + 
  scale_color_manual(values = zeta_colors) +
  ggtitle("AFIO-CAVI") +
  theme_bw()



# Combine plots
p = (p1 + p2) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom") &
  labs(x = "Iteration", y = expression(posterior~of~zeta[t]))

ggsave(plot = p, filename = file.path(save_path, "plots/res_zeta.pdf"),
       device = "pdf", 
       width = 9, height = 6, dpi = 600)

################################################################################
#demonstration of the partial-update idea: increase of ELBO + subsample

#plot the 2 epsilon curves in theorical part
#Full
#Partial Random
#Partial eELBO
#Partial eIter
#Partial eIter adaptive 

#"full", "partial_random","partial_eELBO", "partial_eIter","partial_eIter_adaptive", "full_anneal", "partial_eIter_adaptive_anneal")

df_full = read.table(file.path(save_path, "perform/perform_6.csv"))
df_partial_random = read.table(file.path(save_path, "perform/perform_106.csv"))
df_partial_eELBO = read.table(file.path(save_path, "perform/perform_206.csv"))
df_partial_eIter = read.table(file.path(save_path, "perform/perform_306.csv"))
df_partial_eIter_adaptive = read.table(file.path(save_path, "perform/perform_406.csv"))
df_partial_eIter_adaptive_anneal = read.table(file.path(save_path, "perform/perform_506.csv"))


# Set common theme with shared y-axis limits
all_elbo <- c(
  df_full$ELBO,
  df_partial_random$ELBO,
  df_partial_eELBO$ELBO,
  df_partial_eIter$ELBO,
  df_partial_eIter_adaptive$ELBO,
  df_partial_eIter_adaptive_anneal$ELBO
)
y_limits <- range(all_elbo[is.finite(all_elbo)], na.rm = TRUE)

# Theme to remove y-axis and set smaller title font
theme_noy_smalltitle <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 10)
)

# Title-only theme for the first plot (with y-axis)
theme_y_smalltitle <- theme(
  plot.title = element_text(size = 10)
)

p1 <- df_full %>%
  ggplot(aes(x = iter, y = ELBO, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  # scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  ggtitle("Full") + theme_y_smalltitle

p2 <- df_partial_random %>%
  ggplot(aes(x = iter, y = ELBO, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  # scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  ggtitle("Partial Random") + theme_noy_smalltitle

p3 <- df_partial_eELBO %>%
  ggplot(aes(x = iter, y = ELBO, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  # scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  ggtitle("Partial eELBO") + theme_noy_smalltitle

p4 <- df_partial_eIter %>%
  ggplot(aes(x = iter, y = ELBO, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  # scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  ggtitle("Partial eIter") + theme_noy_smalltitle

p5 <- df_partial_eIter_adaptive %>%
  ggplot(aes(x = iter, y = ELBO, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  # scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  ggtitle("Partial eIter Adapt") + theme_noy_smalltitle


p6 <- df_partial_eIter_adaptive_anneal %>%
  ggplot(aes(x = iter, y = ELBO, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  # scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  ggtitle("Partial eIter Adapt Anneal") + theme_noy_smalltitle


p = (p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom") &
  labs(x = "Iteration", y = "ELBO")

ggsave(plot = p, filename = file.path(save_path, "plots/res_elbo.pdf"),
       device = "pdf", 
       width = 14, height = 4, dpi = 1200)


## Change of subsample size 

# Get shared y-axis range across all plots (filter out non-finite just in case)
all_subsample_sizes <- c(
  df_partial_random$subsample_size,
  df_partial_eELBO$subsample_size,
  df_partial_eIter$subsample_size
)
y_limits <- range(all_subsample_sizes[is.finite(all_subsample_sizes)], na.rm = TRUE)

# Define theme helpers
theme_y_smalltitle <- theme(
  plot.title = element_text(size = 10)
)

theme_noy_smalltitle <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.title = element_text(size = 10)
)

# Define plots
p2 <- df_partial_random %>%
  ggplot(aes(x = iter, y = subsample_size, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  ggtitle("Partial Random") + theme_y_smalltitle

p3 <- df_partial_eELBO %>%
  ggplot(aes(x = iter, y = subsample_size, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  ggtitle("Partial eELBO") + theme_noy_smalltitle

p4 <- df_partial_eIter %>%
  ggplot(aes(x = iter, y = subsample_size, color = partial)) +
  geom_point() + theme_bw() + ylim(y_limits) +
  ggtitle("Partial eIter") + theme_noy_smalltitle

# Combine plots
(p2 + p3 + p4) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position = "bottom") &
  labs(x = "Iteration", y = "Number of responses updated")


################################################################################
#plot the logistic function
logistic_function = function(x) {
  1/ (1 + exp(-x))
}


x_vals <- seq(0.001, 100, 0.001)
df <- data.frame(
  x = x_vals,
  y = logistic_function(log(x_vals))
)

# Logistic function
p1 = ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "blue", size = 1) +
  theme_minimal() +
  labs(
    y = expression(epsilon),
    x = "ELBO difference",
    title = "Logistic Function"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")


#geometric decay
alpha_ls <- c(0.9, 0.95, 0.99)
df <- expand.grid(i = 1:500, a = alpha_ls) %>%
  mutate(y = a^(i - 1))

# Convert `a` to a factor for labeling
df$a <- factor(df$a)

# Plot
p2 = ggplot(df, aes(x = i, y = y, color = a)) +
  geom_line(size = 0.8) +
  theme_minimal() +
  labs(
    title = expression(paste("Geometric Decay: ", epsilon[i] == alpha^{i-1})),
    x = "i",
    y = expression(epsilon[i]),
    color = expression(alpha)
  )
p1+p2

#truncated AUC

