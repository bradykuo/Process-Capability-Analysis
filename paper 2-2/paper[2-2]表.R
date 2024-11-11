library(stats)

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

calculate_lcb <- function(n, Cpk_hat, xi_hat=1.0, gamma=0.95) {
  result <- uniroot(integrate_func, c(0, Cpk_hat), 
                    n=n, Cpk_hat=Cpk_hat, xi_hat=xi_hat, gamma=gamma)
  return(result$root)
}

start_time <- Sys.time()
n_values <- seq(10, 200, by = 5)
c_pk_hat_values_A <- seq(0.7, 1.8, by = 0.1)
c_pk_hat_values_B <- seq(1.9, 3.0, by = 0.1)

# Calculate Panel A
panel_A <- sapply(c_pk_hat_values_A, function(c_pk) {
  sapply(n_values, function(n) round(calculate_lcb(n, c_pk), 3))
})

# Calculate Panel B
panel_B <- sapply(c_pk_hat_values_B, function(c_pk) {
  sapply(n_values, function(n) round(calculate_lcb(n, c_pk), 3))
})

# Create data frames
result_df_A <- data.frame(n = n_values)
result_df_A <- cbind(result_df_A, panel_A)

result_df_B <- data.frame(n = n_values)
result_df_B <- cbind(result_df_B, panel_B)

# Set column names
colnames(result_df_A) <- c("n", as.character(c_pk_hat_values_A))
colnames(result_df_B) <- c("n", as.character(c_pk_hat_values_B))

# Print the tables
cat("Panel A:\n")
print(result_df_A, row.names = FALSE)

cat("\nPanel B:\n")
print(result_df_B, row.names = FALSE)

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")
print(total_time)
