library(nleqslv)

# Define the function to calculate the critical value c0
calculate_c0_nleqslv <- function(C, Cp, n, alpha) {
  # Define the objective function for root finding
  objective <- function(c0) {
    integrand <- function(y) {
      G_term <- pchisq((n - 1) * (3 * Cp * sqrt(n) - y)^2 / (9 * n * c0^2), df = n - 1)
      f_term <- dnorm(y + 3 * (Cp - C) * sqrt(n)) + dnorm(y - 3 * (Cp - C) * sqrt(n))
      G_term * f_term
    }
    
    # Handle potential errors in integration
    integral_result <- tryCatch(
      integrate(integrand, lower = 0, upper = 3 * Cp * sqrt(n))$value,
      error = function(e) { NA }
    )
    
    # Ensure valid output
    if (is.na(integral_result)) {
      stop("Integration failed. Check the function or parameter values.")
    }
    
    integral_result - alpha
  }
  
  # Use nleqslv to solve for c0
  result <- nleqslv(x = 1.5, fn = objective)  # Initial guess can be adjusted
  
  # Check the result
  if (!is.null(result$x)) {
    cat("Calculated c0:", result$x, "\n")
    return(result$x)
  } else {
    stop("nleqslv did not return a valid solution.")
  }
}

# Example parameters
C <- 1.33
n <- 50
alpha <- 0.05
Cp <- if (n < 100) C + 0.33 else C + 0.12

# Calculate c0 using nleqslv
c0 <- calculate_c0_nleqslv(C, Cp, n, alpha)
