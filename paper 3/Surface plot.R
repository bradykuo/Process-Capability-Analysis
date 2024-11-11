# Load required libraries
library(ggplot2)
library(pracma)
library(tidyr)
library(dplyr)
library(metR)
library(gridExtra)
library(plotly)

# Define parameters
alpha <- 0.05  # Producer's risk
beta <- 0.1    # Consumer's risk
C_AQL <- 1.33  # Acceptable Quality Level
C_LTPD <- 1.00 # Lot Tolerance Percent Defective
xi <- 1.0      # Using ξ = 1.0 as per paper
b1 <- 3 * C_AQL + abs(xi)
b2 <- 3 * C_LTPD + abs(xi)

# Generate grid points
n_values <- seq(10, 150, length.out = 50)
c0_values <- seq(0.8, 1.5, length.out = 50)
grid <- expand.grid(n = n_values, c0 = c0_values)

# Functions for S1 and S2 remain the same
calculate_S1 <- function(n, c0) {
  integrand <- function(t) {
    G_term <- pchisq((n - 1) * (b1 * sqrt(n) - t)^2 / (9 * n * c0^2), df = n-1)
    phi_term <- dnorm(t + xi*sqrt(n)) + dnorm(t - xi*sqrt(n))
    return(G_term * phi_term)
  }
  
  upper_limit <- b1 * sqrt(n)
  result <- try(integrate(integrand, lower = 0, upper = upper_limit)$value - (1 - alpha), silent = TRUE)
  
  if(inherits(result, "try-error")) {
    return(NA)
  }
  return(result)
}

calculate_S2 <- function(n, c0) {
  integrand <- function(t) {
    G_term <- pchisq((n - 1) * (b2 * sqrt(n) - t)^2 / (9 * n * c0^2), df = n-1)
    phi_term <- dnorm(t + xi*sqrt(n)) + dnorm(t - xi*sqrt(n))
    return(G_term * phi_term)
  }
  
  upper_limit <- b2 * sqrt(n)
  result <- try(integrate(integrand, lower = 0, upper = upper_limit)$value - beta, silent = TRUE)
  
  if(inherits(result, "try-error")) {
    return(NA)
  }
  return(result)
}


# Calculate S1 and S2 for each point
grid$S1 <- mapply(calculate_S1, grid$n, grid$c0)
grid$S2 <- mapply(calculate_S2, grid$n, grid$c0)

# Reshape data for surface plots - note the transpose for axis swap
S1_matrix <- matrix(grid$S1, nrow = length(n_values), ncol = length(c0_values))
S2_matrix <- matrix(grid$S2, nrow = length(n_values), ncol = length(c0_values))

# Create surface plot for S1 with swapped axes
surface_plot_S1 <- plot_ly(x = ~c0_values, y = ~n_values, z = ~t(S1_matrix)) %>%
  add_surface(
    colorscale = list(
      c(0, 1),
      c("blue", "lightblue")
    )
  ) %>%
  layout(
    title = "Surface Plot of S1(n, c0)",
    scene = list(
      xaxis = list(title = "c0", range = c(0.8, 1.5)),
      yaxis = list(title = "n", range = c(10, 150)),
      zaxis = list(title = "S1"),
      camera = list(
        eye = list(x = 1.5, y = -1.5, z = 1.5)  # Adjusted camera angle
      )
    )
  )

# Create surface plot for S2 with swapped axes
surface_plot_S2 <- plot_ly(x = ~c0_values, y = ~n_values, z = ~t(S2_matrix)) %>%
  add_surface(
    colorscale = list(
      c(0, 1),
      c("red", "pink")
    )
  ) %>%
  layout(
    title = "Surface Plot of S2(n, c0)",
    scene = list(
      xaxis = list(title = "c0", range = c(0.8, 1.5)),
      yaxis = list(title = "n", range = c(10, 150)),
      zaxis = list(title = "S2"),
      camera = list(
        eye = list(x = 1.5, y = -1.5, z = 1.5)  # Adjusted camera angle
      )
    )
  )

# Create combined surface plot with swapped axes
combined_surface_plot <- plot_ly() %>%
  add_surface(
    x = ~c0_values, 
    y = ~n_values, 
    z = ~t(S1_matrix),
    colorscale = list(
      c(0, 1),
      c("blue", "lightblue")
    ),
    opacity = 0.7,
    name = "S1"
  ) %>%
  add_surface(
    x = ~c0_values, 
    y = ~n_values, 
    z = ~t(S2_matrix),
    colorscale = list(
      c(0, 1),
      c("red", "pink")
    ),
    opacity = 0.7,
    name = "S2"
  ) %>%
  layout(
    title = "Surface Plot of S1 and S2",
    scene = list(
      xaxis = list(title = "c0", range = c(0.8, 1.5)),
      yaxis = list(title = "n", range = c(10, 150)),
      zaxis = list(title = "Value"),
      camera = list(
        eye = list(x = 1.5, y = -1.5, z = 1.5)  # Adjusted camera angle
      )
    )
  )

# Add zero plane to help visualize intersection
combined_surface_plot <- combined_surface_plot %>%
  add_trace(
    type = 'surface',
    x = ~c0_values,
    y = ~n_values,
    z = matrix(0, nrow = length(c0_values), ncol = length(n_values)),
    opacity = 0.3,
    colorscale = list(c(0,1), c("gray", "gray")),
    showscale = FALSE,
    name = "Zero Plane"
  )

# Return the plots in a list
list(
  S1_surface = surface_plot_S1,
  S2_surface = surface_plot_S2,
  combined_surface = combined_surface_plot
)