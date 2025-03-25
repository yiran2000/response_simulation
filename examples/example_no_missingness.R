# library(atlasqtl)
library(ggplot2)
library(dplyr)
library(data.table)
seed <- 123; set.seed(seed)

###################
## Simulate data ##
###################

# Example with small problem sizes:
#
n <- 1000; p <- 200; p_act <- 5; q <- 200; q_act <- 10

# Candidate predictors (subject to selection)
#
# Here example with common genetic variants under Hardy-Weinberg equilibrium
#
X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)

# shuffle indices 
shuff_x_ind <- sample(p)
shuff_y_ind <- sample(q)

X <- cbind(X_act, X_inact)[, shuff_x_ind]

# Association pattern and effect sizes
#
pat <- matrix(FALSE, ncol = q, nrow = p)
bool_x <- shuff_x_ind <= p_act
bool_y <- shuff_y_ind <= q_act

pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = q_act)
pat_act[sample(p_act * q_act, floor(p_act * q_act / 5))] <- 1
beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))

pat[bool_x, bool_y] <- pat_act

# Gaussian responses
#
Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta_act), nrow = n)
Y_inact <- matrix(rnorm(n * (q - q_act)), nrow = n)

Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]

########################
## Infer associations ##
########################

# Expectation and variance for the prior number of predictors associated with
# each response
#
# p0 <- c(mean(colSums(pat)), 10)
mu_t <- 1; v_t <- 4 

devtools::load_all()


system.time(res_atlas <- atlasqtl(as.matrix(Y), as.matrix(X),
                                  p0 = c(mu_t, v_t),
                                  user_seed = 1, maxit= 50000,
                                  batch = "y",
                                  thinned_elbo_eval = F,
                                  # anneal = NULL,
                                  tol_loose = 0.01,
                                  tol_tight = 0.01,
                                  burn_in = 50000,
                                  maxit_full = 1,
                                  maxit_subsample = 1000000,
                                  n_partial_update = 500000,
                                  # iter_ladder = c(5, 10, 15, 20, 25, 30, 40, 60, 80, 100),
                                  # e_ladder = c(0.9, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05),
                                  epsilon= c(2, 1.5, 0.25),
                                  partial_elbo = F, #whether we calculate elbo with all responses or only a part
                                  partial_elbo_eval = F, #whether diff_lb = lb_new -lb_old or (lb_new-lb_old)/length(sample_q)
                                  eval_perform = T))

########################
## Algorithm evaluation ##
########################

################################################################################
#check zeta
zeta_df = do.call(rbind, res_atlas$zeta_ls)


zeta_long <- reshape2::melt(zeta_df, variable.name = "Column", value.name = "Value")
colnames(zeta_long) = c("Iteration","Column","Value")

zeta_long = zeta_long %>% mutate(if_sig_protein = if_else(Column %in% which(colSums(pat)>0), 1, 0))


ggplot(zeta_long, aes(x = Iteration, y = Value, group = Column, color = as.factor(if_sig_protein))) +
  geom_line(alpha = 0.5) +  # Line plot
  geom_point(alpha = 0.7, size = 1) +  # Dot plot
  labs(x = "Iteration",
       y = "Zeta") +
  theme_minimal() 

#check select prob

select_prob_df = do.call(rbind, res_atlas$select_prob_ls)
select_prob_long <- reshape2::melt(select_prob_df, variable.name = "Column", value.name = "Value")
colnames(select_prob_long) = c("Iteration","Column","Value")
select_prob_long = select_prob_long %>% mutate(if_sig_protein = if_else(Column %in% which(colSums(pat)>0), 1, 0))

ggplot(select_prob_long, aes(x = Iteration, y = Value, group = Column, color = as.factor(if_sig_protein))) +
  geom_line(alpha = 0.5) +  # Line plot
  geom_point(alpha = 0.7, size = 1) +  # Dot plot
  labs(x = "Iteration",
       y = "select_prob") +
  theme_minimal() 


#check r_vc

r_vc_df = do.call(rbind, res_atlas$r_vc_ls)
r_vc_long <- reshape2::melt(r_vc_df, variable.name = "Column", value.name = "Value")
colnames(r_vc_long) = c("Iteration","Column","Value")
r_vc_long = r_vc_long %>% mutate(if_sig_protein = if_else(Column %in% which(colSums(pat)>0), 1, 0))

ggplot(r_vc_long, aes(x = Iteration, y = Value, group = Column, color = as.factor(if_sig_protein))) +
  geom_line(alpha = 0.1) +  # Line plot
  geom_point(alpha = 0.7, size = 1) +  # Dot plot
  labs(x = "Iteration",
       y = "r_vc") +
  theme_minimal() 



perform_df= res_atlas$perform_df

# Check how the ELBO changes over iterations 
perform_df %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial)) + geom_point()

# Check how the difference of ELBO changes over iterations 
perform_df %>% 
  ggplot(aes(x = iter, y = ELBO_diff, color = partial)) + geom_point()

# Check how e changes over iterations
perform_df %>% 
  ggplot(aes(x = iter, y = e, color = partial)) + geom_point()


# lognormal_cdf <- function(x, mu, sigma, m) {
#   return(m + (1 - m) * pnorm((log(x) - mu) / sigma))
# }
# 
# perform_df %>% mutate(
#    new_e = lognormal_cdf(ELBO_diff, mu =1, sigma = 1, m = 0.1)
# ) %>% ggplot(aes(x = iter, y =  new_e)) + geom_point()
# 
# lognormal_cdf(perform_df$ELBO_diff,mu=2, sigma = 1, m = 0.25)
# 
# perform_df %>% mutate(
#   new_e = lognormal_cdf(iter, mu = 1, sigma = 0.5)
# ) %>% ggplot(aes(x = iter, y =  new_e, color = subsample)) + geom_point()


#

res_atlas$perform_df %>% 
  as.data.frame() %>% 
  ggplot(aes(x = iter, y = ELBO_diff, color = partial))+ 
  geom_point()+
  geom_line()

res_atlas$perform_df %>% 
  as.data.frame() %>% 
  ggplot(aes(x = iter, y = ELBO, color = partial))+ 
  geom_point()+
  geom_line()
 

res_atlas$perform_df %>% 
  as.data.frame() %>% 
  ggplot(aes(x = iter, y = subsample_size, color = partial))+ 
  geom_point()+
  geom_line()

res_atlas$perform_df %>% 
  as.data.frame() %>% 
  ggplot(aes(x = ELBO_diff, y = subsample_size, color = partial))+ 
  geom_point()+
  geom_line()
