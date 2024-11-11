# Function to calculate PPM bounds for given Cpk
calculate_ppm_bounds <- function(cpk) {
  # For one-sided specification:
  upper_ppm <- (1 - pnorm(3 * cpk)) * 1e6
  
  # For two-sided specification:
  # Lower bound is when process is centered (symmetrical case)
  lower_ppm <- 2 * (1 - pnorm(3 * cpk)) * 1e6
  
  return(c(lower_ppm, upper_ppm))
}

# Define Cpk values
cpk_values <- c(0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.24, 1.25, 
                1.30, 1.33, 1.40, 1.45, 1.50, 1.60, 1.67, 1.70, 1.80, 1.90, 2.00)

# Calculate PPM bounds for each Cpk value
results <- t(sapply(cpk_values, calculate_ppm_bounds))

# Create data frame
cpk_table <- data.frame(
  Cpk = cpk_values,
  Lower_Bound = round(results[,1], 3),
  Upper_Bound = round(results[,2], 3)
)

# Print the results
print(cpk_table, row.names = FALSE)

# Verify some values against the paper
cat("\nVerification of selected values:\n")
cat("For Cpk = 1.33: Expected ≈ 66/33 PPM, Calculated:", 
    round(calculate_ppm_bounds(1.33), 3), "PPM\n")
cat("For Cpk = 1.67: Expected ≈ 0.544/0.272 PPM, Calculated:", 
    round(calculate_ppm_bounds(1.67), 3), "PPM\n")