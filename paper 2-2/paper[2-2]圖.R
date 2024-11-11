library(stats)
library(ggplot2)
library(gridExtra)
library(grid)

integrate_func <- function(C, n, Cpk_hat, xi_hat, gamma) {
  b <- 3*C + abs(xi_hat)
  
  integrand <- function(t) {
    G_val <- pchisq((n-1)*(b*sqrt(n)-t)^2 / (9*n*Cpk_hat^2), df=n-1)
    phi_val <- dnorm(t + xi_hat*sqrt(n)) + dnorm(t - xi_hat*sqrt(n))
    return(G_val * phi_val)
  }
  
  result <- integrate(integrand, lower=0, upper=b*sqrt(n))
  return(result$value - (1 - gamma))
}

calculate_lcb <- function(n, Cpk_hat, xi_hat, gamma=0.95) {
  result <- uniroot(integrate_func, c(0, Cpk_hat), 
                    n=n, Cpk_hat=Cpk_hat, xi_hat=xi_hat, gamma=gamma)
  return(result$root)
}

calculate_C_vs_xi <- function(Cpk_hat, n_values, xi_values) {
  results <- expand.grid(n = n_values, xi = xi_values)
  results$C <- mapply(calculate_lcb, results$n, MoreArgs = list(Cpk_hat = Cpk_hat, gamma = 0.95), xi_hat = results$xi)
  return(results)
}

start_time <- Sys.time()
n_values <- c(30, 50, 70, 100, 150, 200)
xi_values <- seq(0, 3, by = 0.1)
Cpk_hat_values <- c(0.7, 0.9, 1.2, 2.0, 2.5, 3.0)

plots <- list()

for (i in 1:length(Cpk_hat_values)) {
  Cpk_hat <- Cpk_hat_values[i]
  results <- calculate_C_vs_xi(Cpk_hat, n_values, xi_values)
  
  p <- ggplot(results, aes(x = xi, y = C, color = factor(n))) +
    geom_line() +
    labs(title = paste("Ĉpk =", Cpk_hat),
         x = "|ξ|", y = "C") +
    theme_minimal() +
    theme(legend.position = "none") +  # Remove individual legends
    scale_color_discrete(name = "n")
  
  plots[[i]] <- p
}

# Create a common legend
legend_plot <- ggplot(results, aes(x = xi, y = C, color = factor(n))) +
  geom_line() +
  scale_color_discrete(name = "n") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Extract the legend
tmp <- ggplot_gtable(ggplot_build(legend_plot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

# Arrange plots in a grid
arranged_plots <- gridExtra::grid.arrange(
  grobs = c(plots, list(legend)),
  layout_matrix = rbind(matrix(1:6, nrow = 2, byrow = TRUE),
                        c(NA, 7, NA)),
  top = grid::textGrob("C vs |ξ| for different Ĉpk values", gp = grid::gpar(fontsize = 16, font = 2))
)

# Save the plot
ggsave("Cpk_plots.png", arranged_plots, width = 15, height = 10, dpi = 300)

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")
print(total_time)