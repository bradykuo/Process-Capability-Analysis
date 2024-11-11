# Core S1 and S2 functions remain the same
calculate_S1 <- function(n, c0, c_aql, alpha, xi = 1) {
  b1 <- 3 * c_aql + abs(xi)
  
  integrand <- function(t) {
    G_val <- pchisq((n - 1) * (b1 * sqrt(n) - t)^2 / (9 * n * c0^2), df = n - 1)
    phi_val <- dnorm(t + sqrt(n)) + dnorm(t - sqrt(n))
    return(G_val * phi_val)
  }
  
  result <- integrate(integrand, lower = 0, upper = b1 * sqrt(n))$value - (1 - alpha)
  return(result)
}

calculate_S2 <- function(n, c0, c_ltpd, beta, xi = 1) {
  b2 <- 3 * c_ltpd + abs(xi)
  
  integrand <- function(t) {
    G_val <- pchisq((n - 1) * (b2 * sqrt(n) - t)^2 / (9 * n * c0^2), df = n - 1)
    phi_val <- dnorm(t + sqrt(n)) + dnorm(t - sqrt(n))
    return(G_val * phi_val)
  }
  
  result <- integrate(integrand, lower = 0, upper = b2 * sqrt(n))$value - beta
  return(result)
}

# Modified function to find c0 for a given n
find_c0_for_n <- function(n, alpha, beta, c_aql, c_ltpd) {
  objective <- function(c0) {
    s1 <- calculate_S1(n, c0, c_aql, alpha)
    s2 <- calculate_S2(n, c0, c_ltpd, beta)
    return(abs(s1) + abs(s2))  # Minimize total error
  }
  
  # Use optimize instead of uniroot to minimize total error
  result <- optimize(objective, interval = c(c_ltpd, c_aql + 0.5))
  return(result$minimum)
}

# Function to evaluate error for a given n
evaluate_n <- function(n, alpha, beta, c_aql, c_ltpd) {
  c0 <- find_c0_for_n(n, alpha, beta, c_aql, c_ltpd)
  s1 <- calculate_S1(n, c0, c_aql, alpha)
  s2 <- calculate_S2(n, c0, c_ltpd, beta)
  return(abs(s1) + abs(s2))
}

# Main function to find n and c0
find_n_c0 <- function(alpha, beta, c_aql, c_ltpd) {
  # Set search range 
  n_min <- 30
  n_max <- 1000
  
  best_n <- NA
  best_c0 <- NA
  min_error <- Inf
  
  for(n in n_min:n_max) {
    c0 <- find_c0_for_n(n, alpha, beta, c_aql, c_ltpd)
    error <- evaluate_n(n, alpha, beta, c_aql, c_ltpd)
    
    if(error < min_error) {
      min_error <- error
      best_n <- n
      best_c0 <- c0
    }
  }
  
  return(c(n = best_n, c0 = round(best_c0, 4)))
}

# Function to generate Table 3
generate_table3 <- function() {
  # Define parameters
  alpha_values <- c(0.01, 0.025, 0.05, 0.075, 0.10)
  beta_values <- c(0.01, 0.025, 0.05, 0.075, 0.10)
  cases <- list(
    c(c_aql = 1.33, c_ltpd = 1.00),
    c(c_aql = 1.50, c_ltpd = 1.33),
    c(c_aql = 1.67, c_ltpd = 1.33),
    c(c_aql = 2.00, c_ltpd = 1.67)
  )
  
  # Create results matrix
  results <- matrix(NA, nrow = length(alpha_values) * length(beta_values), 
                    ncol = 8 + 2)  # alpha, beta + 4 pairs of (n, c0)
  colnames(results) <- c("alpha", "beta", 
                         "n1", "c01", "n2", "c02", "n3", "c03", "n4", "c04")
  
  row <- 1
  for(alpha in alpha_values) {
    for(beta in beta_values) {
      results[row, 1:2] <- c(alpha, beta)
      
      for(i in 1:length(cases)) {
        case <- cases[[i]]
        result <- find_n_c0(alpha, beta, case["c_aql"], case["c_ltpd"])
        results[row, (2*i+1):(2*i+2)] <- result
      }
      
      row <- row + 1
      cat(sprintf("Calculated: alpha = %.3f, beta = %.3f\n", alpha, beta))
    }
  }
  
  # Convert to data frame
  results_df <- as.data.frame(results)
  
  # Print formatted table
  cat("\nTable 3: Required sample sizes (n) and critical acceptance values (c0)\n")
  cat("α      β       CAQL=1.33      CAQL=1.50      CAQL=1.67      CAQL=2.00\n")
  cat("               CLTPD=1.00      CLTPD=1.33      CLTPD=1.33      CLTPD=1.67\n")
  cat("               n      c0       n      c0       n      c0       n      c0\n")
  cat(rep("-", 75), "\n")
  
  for(i in 1:nrow(results_df)) {
    cat(sprintf("%.3f  %.3f  ", results_df[i,1], results_df[i,2]))
    for(j in 1:4) {
      cat(sprintf("%3d  %.4f  ", results_df[i,2*j+1], results_df[i,2*j+2]))
    }
    cat("\n")
  }
  
  return(results_df)
}

# Generate Table 3
start_time <- proc.time()
table3 <- generate_table3()
end_time <- proc.time()
run_time <- end_time - start_time
print(run_time)