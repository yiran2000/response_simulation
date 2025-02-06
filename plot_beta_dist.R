# Function to plot Beta distribution
plot_beta <- function(alpha, beta, x_min = 0, x_max = 1, n_points = 1000) {
  # Input validation
  if (alpha <= 0 || beta <= 0) {
    stop("Alpha and Beta must be positive values.")
  }
  
  # Generate x values
  x <- seq(x_min, x_max, length.out = n_points)
  
  # Compute Beta density
  y <- dbeta(x, shape1 = alpha, shape2 = beta)
  
  # Create the plot
  plot(x, y, type = "l", col = "blue", lwd = 2,
       main = paste("Beta Distribution (alpha =", alpha, ", beta =", beta, ")"),
       xlab = "x", ylab = "Density",
       ylim = c(0, max(y) * 1.1))
  
  # Add a grid for better visualization
  grid()
}

# Example usage:
# Plot Beta(2, 5)
plot_beta(2, 100)
