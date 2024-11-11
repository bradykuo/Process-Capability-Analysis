library(parallel)
library(knitr)
library(kableExtra)
library(Rcpp)

# Rcpp implementation of Gpk calculation
cppFunction('
NumericVector calculate_Gpk_cpp(double x_bar, double s2, int n, double L, double U, int num_simulations) {
  double d = (U - L) / 2;
  double M = (U + L) / 2;
  
  NumericVector Z = rnorm(num_simulations);
  NumericVector U2 = rchisq(num_simulations, n-1);
  
  NumericVector T_mu = x_bar - sqrt((n-1.0)/n) * (Z * sqrt(s2) / sqrt(n));
  NumericVector T_sigma2 = s2 * (n-1) / U2;
  
  return (d - abs(T_mu - M)) / (3 * sqrt(T_sigma2));
}
')

calculate_Cpk <- function(mu, sigma, L, U) {
  min((U - mu) / (3 * sigma), (mu - L) / (3 * sigma))
}

calculate_other_limits <- function(x_bar, s, n, L, U, alpha) {
  Cpk_hat <- calculate_Cpk(x_bar, s, L, U)
  z <- qnorm(1 - alpha)
  
  Bpk <- Cpk_hat - z * sqrt(1 / (9 * n) + Cpk_hat^2 / (2 * (n - 1)))
  KHpk <- Cpk_hat * (1 - z / sqrt(2 * (n - 1)))
  Hpk <- Cpk_hat - z * sqrt((n - 1) / (9 * n * (n - 3)) + Cpk_hat^2 * (1 / (2 * (n - 3)) * (1 + 6 / (n - 1))))
  Npk <- sqrt(1 - 2 / (5 * (n - 1))) * Cpk_hat - z * sqrt(Cpk_hat^2 / (2 * (n - 1)) + 1 / (9 * n))
  
  c(Bpk, KHpk, Hpk, Npk)
}

simulate_Cpk <- function(n, mu, sigma, L, U, Cpk, alpha, num_outer_simulations, num_inner_simulations) {
  true_Cpk <- Cpk
  
  x <- matrix(rnorm(n * num_outer_simulations, mu, sigma), ncol = n)
  x_bar <- rowMeans(x)
  s <- apply(x, 1, sd)
  s2 <- s^2
  
  Gpk_values <- sapply(seq_len(num_outer_simulations), function(i) {
    quantile(calculate_Gpk_cpp(x_bar[i], s2[i], n, L, U, num_inner_simulations), alpha)
  })
  
  other_limits <- t(sapply(seq_len(num_outer_simulations), function(i) {
    calculate_other_limits(x_bar[i], s[i], n, L, U, alpha)
  }))
  
  limits <- cbind(Gpk_values, other_limits)
  
  coverage <- colMeans(limits <= true_Cpk)
  expected_values <- colMeans(limits)
  
  list(coverage = round(coverage, 4), expected_values = round(expected_values, 4))
}

format_results <- function(results) {
  formatted_output <- do.call(rbind, lapply(names(results), function(Cpk) {
    do.call(rbind, lapply(names(results[[Cpk]]), function(n) {
      do.call(rbind, lapply(names(results[[Cpk]][[n]]), function(alpha) {
        c(Cpk, n, 1-as.numeric(alpha),
          results[[Cpk]][[n]][[alpha]]$coverage,
          results[[Cpk]][[n]][[alpha]]$expected_values)
      }))
    }))
  }))
  
  formatted_output <- as.data.frame(formatted_output)
  colnames(formatted_output) <- c("Cpk", "n", "1-α", 
                                  "Gpk", "Bpk", "KHpk", "Hpk", "Npk",
                                  "E(Gpk)", "E(Bpk)", "E(KHpk)", "E(Hpk)", "E(Npk)")
  
  # Convert numeric columns
  formatted_output$Cpk <- as.numeric(as.character(formatted_output$Cpk))
  formatted_output$n <- as.numeric(as.character(formatted_output$n))
  formatted_output$`1-α` <- as.numeric(as.character(formatted_output$`1-α`))
  
  numeric_cols <- 4:13  # columns 4 through 13 should be numeric
  formatted_output[numeric_cols] <- lapply(formatted_output[numeric_cols], 
                                           function(x) round(as.numeric(as.character(x)), 4))
  
  return(formatted_output)
}

# Main simulation
set.seed(123)  # for reproducibility
Cpk_values <- c(1, 1.33, 1.5, 2, 2.5, 3)
n_values <- c(10, 20, 30, 40, 50)
alpha_values <- c(0.1, 0.05)
num_outer_simulations <- 10000
num_inner_simulations <- 10000
L <- 7
U <- 14
mu <- 10
d <- (U - L) / 2
M <- (U + L) / 2

start_time <- Sys.time()

results <- mclapply(Cpk_values, function(Cpk) {
  sigma <- (d - abs(mu - M)) / (3 * Cpk)
  
  lapply(n_values, function(n) {
    lapply(alpha_values, function(alpha) {
      cat(sprintf("Processing Cpk = %.2f, n = %d, alpha = %.2f\n", Cpk, n, alpha))
      simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, num_outer_simulations, num_inner_simulations)
    })
  })
}, mc.cores = detectCores())

names(results) <- as.character(Cpk_values)
for (i in seq_along(results)) {
  names(results[[i]]) <- as.character(n_values)
  for (j in seq_along(results[[i]])) {
    names(results[[i]][[j]]) <- as.character(alpha_values)
  }
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")

formatted_results <- format_results(results)

# Modify column names to add spacing
names(formatted_results)[3] <- "1-α   "
names(formatted_results)[8] <- "Npk   "

print(kable(formatted_results, 
            format = "pipe",
            digits = 4,
            align = c('r', 'r', 'r', rep('r', 5), rep('r', 5))) %>%
        add_header_above(c(" " = 3, 
                           "Coverage probability" = 5,
                           "Expected value" = 5)) %>%
        add_header_above(c(" " = 13)) %>%
        kable_styling(full_width = F))

# Print execution time
cat(sprintf("\nTotal execution time: %.2f seconds\n", total_time))