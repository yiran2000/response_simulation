#plot beta distribution

x = seq(0, 1, 0.001)
p3 = data.frame(
  x = x,
  y = c(dbeta(x, shape1 = 1, shape2 = 10),dbeta(x, shape1 = 1, shape2 = 40)),
  rho = as.factor(c(rep(10, length(x)), rep(40, length(x))))
) %>% ggplot(aes(x = x, y = y, color = rho))+
  geom_point(size = 0.5)+
  theme_bw()+ggtitle("Density of Beta(1, rho)")

p4 = data.frame(
  x = x,
  y = c(qbeta(x, shape1 = 1, shape2 = 10),qbeta(x, shape1 = 1, shape2 = 40)),
  rho = as.factor(c(rep(10, length(x)), rep(40, length(x))))
) %>% ggplot(aes(x = x, y = y, color = rho))+geom_point(size = 0.5)+
  theme_bw()+ggtitle("Distribution of Beta(1, rho)")
p5 = p3+p4+plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(filename = "Beta_distribution.pdf", plot = p5, width = 7, height = 4)