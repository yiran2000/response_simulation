# Load ggplot2 for visualization
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

emax <- function(x, e0=0, emax = 1, ec50, n) {
  e0 + ((emax-e0) * x^n) / (ec50^n +x^n)
}

x =seq(0, 100, 0.1)

# Define parameter combinations
param_grid <- expand.grid(ec50 = 5, n = c(1, 2, 5))

# Compute y values for each parameter set
df <- param_grid %>%
  mutate(y_values = map2(ec50, n, ~emax(x, e0 = 0.2, ec50 = .x, n = .y))) %>%
  unnest(y_values) %>%
  mutate(x = rep(x, times = nrow(param_grid)))


ggplot(df %>% filter(x<50), aes(x = x, y = y_values, color = factor(n))) +
  geom_line() +
  theme_minimal()



