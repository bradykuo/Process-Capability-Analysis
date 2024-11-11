# Load required libraries
library(ggplot2)
library(dplyr)
library(pracma)

# Function to calculate the integral in equations (9) and (10)
calculate_integral <- function(n, c0, C_value, xi) {
  b <- 3 * C_value + abs(xi)
  
  integrand <- function(t) {
    # G is the chi-square CDF with n-1 degrees of freedom
    G_term <- pchisq((n - 1) * (b * sqrt(n) - t)^2 / (9 * n * c0^2), df = n-1)
    # phi is the standard normal PDF
    phi_term <- dnorm(t + xi*sqrt(n)) + dnorm(t - xi*sqrt(n))
    return(G_term * phi_term)
  }
  
  # Numerical integration
  result <- integrate(integrand, lower = 0, upper = b*sqrt(n))$value
  return(result)
}

# Function to find c0 for given n and parameters
find_c0 <- function(n, C_AQL, C_LTPD, alpha, beta, xi) {
  objective <- function(c0) {
    eq1 <- calculate_integral(n, c0, C_AQL, xi) - (1 - alpha)
    eq2 <- calculate_integral(n, c0, C_LTPD, xi) - beta
    return(eq1^2 + eq2^2)
  }
  
  # Find optimal c0
  result <- optimize(objective, interval = c(0.8, 2.0))
  return(result$minimum)
}

# Generate data for plots
xi_values <- seq(0, 2, by = 0.1)
CAQL_values <- c(1.33, 1.50, 1.67, 2.00)
CLTPD <- 1.00
alpha <- 0.05
beta <- 0.05

# Initialize data frames
sample_size_data <- data.frame()
critical_value_data <- data.frame()

# Calculate values for each CAQL and xi combination
for(CAQL in CAQL_values) {
  for(xi in xi_values) {
    # Find minimum n that satisfies both equations
    n_min <- 10  # Start with small n
    while(TRUE) {
      c0 <- find_c0(n_min, CAQL, CLTPD, alpha, beta, xi)
      
      # Check if equations are satisfied
      eq1 <- calculate_integral(n_min, c0, CAQL, xi) >= (1 - alpha)
      eq2 <- calculate_integral(n_min, c0, CLTPD, xi) <= beta
      
      if(eq1 && eq2) break
      n_min <- n_min + 1
      if(n_min > 100) break  # Safety break
    }
    
    # Store results
    sample_size_data <- rbind(sample_size_data, 
                              data.frame(xi = xi, n = n_min, CAQL = as.factor(CAQL)))
    critical_value_data <- rbind(critical_value_data, 
                                 data.frame(xi = xi, c0 = c0, CAQL = as.factor(CAQL)))
  }
}

# Create plot (a)
plot_a <- ggplot(sample_size_data, aes(x = xi, y = n, color = CAQL)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Required Sample Size vs ξ",
    x = "ξ",
    y = "Required sample size n"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80"),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# Create plot (b)
plot_b <- ggplot(critical_value_data, aes(x = xi, y = c0, color = CAQL)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(1.0, 1.6)) +
  labs(
    title = "Critical Acceptance Value vs ξ",
    x = "ξ",
    y = "Critical value Co"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80"),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# Print plots
print(plot_a)
print(plot_b)

# Save plots if needed
# ggsave("sample_size_plot.png", plot_a, width = 8, height = 6)
# ggsave("critical_value_plot.png", plot_b, width = 8, height = 6)