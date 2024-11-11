# Function to calculate PPM bounds for given Cpk
calculate_ppm_bounds <- function(cpk) {
  # For one-sided specification:
  upper_ppm <- (1 - pnorm(3 * cpk)) * 1e6
  
  # For two-sided specification:
  # Lower bound is when process is centered (symmetrical case)
  lower_ppm <- 2 * (1 - pnorm(3 * cpk)) * 1e6
  
  return(c(lower_ppm, upper_ppm))
}

# Create sequence of Cpk values from 0.6 to 2.0
cpk_seq <- seq(0.6, 2.0, by = 0.01)

# Calculate PPM bounds for each Cpk value
results <- t(sapply(cpk_seq, calculate_ppm_bounds))

# Set up the plotting parameters
par(mar = c(5, 5, 2, 2))  # Adjust margins

# Create the plot
plot(cpk_seq, results[,1], type = "l", 
     xlab = expression(C[pk]), 
     ylab = "PPM",
     ylim = c(0, 80000),        # Set y-axis limits exactly as in paper
     xlim = c(0.6, 2.0),        # Set x-axis limits exactly as in paper
     yaxt = "n",                # Remove default y-axis
     cex.lab = 1.2,            # Adjust label size
     lwd = 1)                  # Thin lines as in paper

# Add upper bound line
lines(cpk_seq, results[,2], lwd = 1)

# Create custom y-axis with specific breaks
axis(2, at = seq(0, 80000, by = 20000), las = 1)

# Add title below the plot
title(sub = expression(paste("Fig. 1. The bounds on nonconforming units in PPM versus ", C[pk], ".")),
      line = 3.5)