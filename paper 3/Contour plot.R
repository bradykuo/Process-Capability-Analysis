# Load required libraries
library(ggplot2)
library(pracma)
library(tidyr)
library(dplyr)
library(metR)
library(gridExtra)

# Define parameters
alpha <- 0.05  # Producer's risk
beta <- 0.1    # Consumer's risk
C_AQL <- 1.33  # Acceptable Quality Level
C_LTPD <- 1.00 # Lot Tolerance Percent Defective
xi <- 1.0      # Using ξ = 1.0 as per paper
b1 <- 3 * C_AQL + abs(xi)
b2 <- 3 * C_LTPD + abs(xi)

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

# Function to find c0 value where S1 = 0 for a given n
find_c0_for_S1 <- function(n) {
  f <- function(c0) calculate_S1(n, c0)
  tryCatch({
    result <- uniroot(f, interval = c(0.8, 1.5), tol = 1e-6)
    return(result$root)
  }, error = function(e) NA)
}

# Function to find c0 value where S2 = 0 for a given n
find_c0_for_S2 <- function(n) {
  f <- function(c0) calculate_S2(n, c0)
  tryCatch({
    result <- uniroot(f, interval = c(0.8, 1.5), tol = 1e-6)
    return(result$root)
  }, error = function(e) NA)
}

# Function to find where S1 and S2 curves intersect
find_intersection_precise <- function() {
  # Function to minimize: difference between c0 values where S1=0 and S2=0
  f <- function(n) {
    n <- ceiling(n)  # Ensure n is integer with ceiling
    c0_s1 <- find_c0_for_S1(n)
    c0_s2 <- find_c0_for_S2(n)
    if(is.na(c0_s1) || is.na(c0_s2)) return(NA)
    return(c0_s1 - c0_s2)
  }
  
  # Find intersection using uniroot
  result <- tryCatch({
    n_intersection <- uniroot(f, interval = c(10, 150), tol = 1e-6)
    n <- ceiling(n_intersection$root)  # Ensure final n is integer with ceiling
    c0 <- find_c0_for_S1(n)  # Recalculate c0 with ceiled n
    return(c(n = n, c0 = c0))
  }, error = function(e) NULL)
  
  return(result)
}

# Find precise intersection point
intersection <- find_intersection_precise()

# Generate grid points with better resolution around the intersection
n_values <- seq(10, 150, length.out = 75)  # Increased resolution
c0_values <- seq(0.8, 1.5, length.out = 75)  # Increased resolution
grid <- expand.grid(n = n_values, c0 = c0_values)

# Calculate S1 and S2 for each point
grid$S1 <- mapply(calculate_S1, grid$n, grid$c0)
grid$S2 <- mapply(calculate_S2, grid$n, grid$c0)

# Create improved S1 contour plot
contour_plot_S1 <- ggplot(grid, aes(x = n, y = c0, z = S1)) +
  geom_contour(aes(z = S1), 
               bins = 10,  # Increased number of contours
               color = "blue") +
  geom_contour(aes(z = S1), 
               breaks = 0,  # Add zero contour line
               color = "blue",
               size = 1) +
  geom_text_contour(aes(z = S1), 
                    skip = 2,  # Adjusted skip for better label spacing
                    rotate = TRUE,
                    size = 3,
                    stroke = 0.2) +
  scale_x_continuous(limits = c(10, 150), 
                     breaks = seq(0, 150, by = 25)) +
  scale_y_continuous(limits = c(0.8, 1.5), 
                     breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "Contour Plot of S1(n, c0)",
    x = "n",
    y = "c0"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Create improved S2 contour plot
contour_plot_S2 <- ggplot(grid, aes(x = n, y = c0, z = S2)) +
  geom_contour(aes(z = S2), 
               bins = 10,  # Increased number of contours
               color = "red") +
  geom_contour(aes(z = S2), 
               breaks = 0,  # Add zero contour line
               color = "red",
               size = 1) +
  geom_text_contour(aes(z = S2), 
                    skip = 2,  # Adjusted skip for better label spacing
                    rotate = TRUE,
                    size = 3,
                    stroke = 0.2) +
  scale_x_continuous(limits = c(10, 150), 
                     breaks = seq(0, 150, by = 25)) +
  scale_y_continuous(limits = c(0.8, 1.5), 
                     breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "Contour Plot of S2(n, c0)",
    x = "n",
    y = "c0"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Function to find intersection point
find_intersection <- function(grid) {
  zero_threshold <- 0.01
  intersection_points <- grid[abs(grid$S1) < zero_threshold & abs(grid$S2) < zero_threshold, ]
  
  if(nrow(intersection_points) > 0) {
    closest_point <- intersection_points[which.min(intersection_points$S1^2 + intersection_points$S2^2), ]
    # Ensure n is integer with ceiling
    closest_point$n <- ceiling(closest_point$n)
    # Recalculate c0 with ceiled n
    c0_value <- find_c0_for_S1(closest_point$n)
    return(c(n = closest_point$n, c0 = c0_value))
  } else {
    return(NULL)
  }
}

# Find intersection point
intersection <- find_intersection(grid)

# Create combined contour plot with intersection point
combined_plot <- ggplot(grid) +
  # Add S1 contours in blue
  geom_contour(aes(x = n, y = c0, z = S1), 
               color = "blue",
               linetype = "solid") +
  geom_text_contour(aes(x = n, y = c0, z = S1),
                    skip = 1,
                    rotate = TRUE,
                    size = 3,
                    color = "blue",
                    stroke = 0.2) +
  
  # Add S2 contours in red
  geom_contour(aes(x = n, y = c0, z = S2),
               color = "red",
               linetype = "dashed") +
  geom_text_contour(aes(x = n, y = c0, z = S2),
                    skip = 1,
                    rotate = TRUE,
                    size = 3,
                    color = "red",
                    stroke = 0.2) +
  
  # Add lines for zero contours to highlight intersection
  geom_contour(aes(x = n, y = c0, z = S1),
               breaks = 0,
               color = "blue",
               size = 1) +
  geom_contour(aes(x = n, y = c0, z = S2),
               breaks = 0,
               color = "red",
               size = 1) +
  
  # Add intersection point and label
  geom_point(data = data.frame(n = intersection["n"], c0 = intersection["c0"]),
             aes(x = n, y = c0),
             color = "black",
             size = 3) +
  annotate("text", 
           x = intersection["n"], 
           y = intersection["c0"] + 0.05,
           label = sprintf("(n, c0) = (%.2f, %.4f)", 
                           intersection["n"], 
                           intersection["c0"]),
           size = 3.5,
           hjust = 0) +
  
  # Customize scales and labels
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 50)) +
  scale_y_continuous(limits = c(0.8, 1.5), breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "Combined Contour Plot of S1(n, c0) and S2(n, c0)",
    subtitle = "S1: Blue (solid), S2: Red (dashed)",
    x = "n",
    y = "C0"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Print intersection point values
print(sprintf("Intersection point: n = %d, c0 = %.4f", 
              round(intersection[1]), intersection[2]))

# Arrange all plots in a grid
grid.arrange(
  contour_plot_S1, 
  contour_plot_S2, 
  combined_plot,
  layout_matrix = rbind(c(1,2), c(3,3)),
  heights = c(1, 1)
)

# Save the combined figure if needed
# ggsave("all_contour_plots.png", arrangeGrob(contour_plot_S1, contour_plot_S2, combined_plot, 
#        layout_matrix = rbind(c(1,2), c(3,3))), 
#        width = 12, height = 10, dpi = 300)