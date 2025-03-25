library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(patchwork)
library(stringr)


simulation_tool_path = "~/project/response_simulation"
data_path = "/rds/user/yl2021/hpc-work/hotspot_sim"
save_path = "/rds/user/yl2021/hpc-work/hotspot_sim/simulation_res"

# simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
# data_path = "data"


source("~/project/response_simulation/data_simulation_tools/simulate_withRealGenome.R")
source("~/project/response_simulation/reshape_atlasQTL_res.R")

qt_ls = c(10, 50, 500, 1000) #set different number of active proteins
m_ls = c(0.01, 0.005, 0.001)

param_table = expand.grid(qt = qt_ls, m = m_ls)

files = lapply(gtools::mixedsort(list.files(path = save_path, pattern = "\\.csv$", full.names = TRUE)), fread)
nc = lapply(files, ncol) %>% unlist()

output = files[which(nc == 21)] %>% rbindlist()


#calculate mean and CI
plot_metric_CI = function(output, metric_col){
  
  metric = sym(metric_col)  # Convert string to symbol
  
  output_summary = output %>%
    group_by(m, qt) %>%
    summarise(
      mean_value = mean(!!metric, na.rm = TRUE),
      sd_value = sd(!!metric, na.rm = TRUE),
      n = n(),
      se = sd_value / sqrt(n),  # Standard Error
      ci_lower = mean_value - qt(0.975, df = n-1) * se,
      ci_upper = mean_value + qt(0.975, df = n-1) * se,
      .groups = "drop"  # Suppress summarise() warning
    )
  
  # Fix aes() by using !!metric instead of hardcoded column
  p = ggplot(output_summary, aes(x = qt, y = mean_value, color = as.factor(m))) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    theme_bw()+
    # facet_wrap(~m)+
    scale_x_continuous(breaks = unique(output_summary$qt))+
    # ggtitle(as.character(metric_col)) +
    labs(color = "m",y = metric_col)
  
  return(list(output_summary, p))
}

p1 = plot_metric_CI(output %>% filter(FP_protein<100), "power_protein")[[2]]
p2 = plot_metric_CI(output %>% filter(FP_protein<100),"FDR_protein")[[2]]
p2 = plot_metric_CI(output %>% filter(FP_protein<100),"FPR_protein")[[2]]

p3 = plot_metric_CI(output %>% filter(FP_protein<100),"FP_protein")[[2]]
p4 = plot_metric_CI(output %>% filter(FP_protein<100),"TP_protein")[[2]]



p = ((p1+labs(x = "Number of active responses", y = "Power")) + 
          (p2+labs(x = "Number of active responses", y = "False discovery rate")) +
          (p3+labs(x = "Number of active responses", y = "False positive rate"))) + 
  plot_layout(nrow = 1, guides = "collect") +  # Collect legend
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave(filename = file.path(save_path, "plots/simulation_comb.pdf"), plot = p, device = "pdf", 
       width = 9, height = 5, dpi = 1200)

