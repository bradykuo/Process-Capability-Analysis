library(rootSolve)

# Define the function to calculate the critical value c0
calculate_c0_multiroot <- function(C, Cp, n, alpha) {
  # Define the objective function for the root finding
  objective <- function(c0) {
    integrand <- function(y) {
      G_term <- pchisq((n - 1) * (3 * Cp * sqrt(n) - y)^2 / (9 * n * c0^2), df = n - 1)
      f_term <- dnorm(y + 3 * (Cp - C) * sqrt(n)) + dnorm(y - 3 * (Cp - C) * sqrt(n))
      G_term * f_term
    }
    integral <- integrate(integrand, lower = 0, upper = 3 * Cp * sqrt(n))$value
    integral - alpha
  }
  
  # Use multiroot to solve for c0
  result <- multiroot(f = objective, start = 1.5)  # Initial guess can be adjusted
  result$root
}

# Example parameters
C <- 1.33
n <- 50
alpha <- 0.05
Cp <- if (n < 100) C + 0.33 else C + 0.12

# Calculate c0 using multiroot
c0 <- calculate_c0_multiroot(C, Cp, n, alpha)
print(c0)
