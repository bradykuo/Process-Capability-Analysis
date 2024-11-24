# Define constants
USL <- 61
LSL <- 40
T <- 49
alpha <- 0.05

mu_list <- c(50, 52)
sigma_list <- c(2.0, 3.0, 3.7)
n_list <- c(20, 40, 70)

B <- 1000  # Bootstrap resamples
N <- 1000  # Simulation replications
z_value <- qnorm(1 - alpha)

# Main simulation loop
for(mu_idx in 1:length(mu_list)) {
  mu <- mu_list[mu_idx]
  
  for(sigma_idx in 1:length(sigma_list)) {
    sigma <- sigma_list[sigma_idx]
    
    # Calculate true values for this combination
    true_cp <- (USL - LSL)/(6 * sigma)
    true_cpk <- min((USL - mu)/(3 * sigma), (mu - LSL)/(3 * sigma))
    true_cpm <- (USL - LSL)/(6 * sqrt(sigma^2 + (mu - T)^2))
    
    for(n_idx in 1:length(n_list)) {
      n <- n_list[n_idx]
      
      # Initialize counters
      sb_counts <- numeric(3)  # [cp, cpk, cpm]
      pb_counts <- numeric(3)
      bcpb_counts <- numeric(3)
      
      for(sim in 1:N) {
        # Generate original sample
        sample_fixed <- rnorm(n, mean = mu, sd = sigma)
        
        # Calculate original statistics
        x_bar_original <- mean(sample_fixed)
        s_original <- sd(sample_fixed)
        temp_original <- sqrt(s_original^2 + (x_bar_original - T)^2)
        
        cp_original <- (USL - LSL)/(6 * s_original)
        cpk_original <- min((USL - x_bar_original)/(3 * s_original),
                            (x_bar_original - LSL)/(3 * s_original))
        cpm_original <- (USL - LSL)/(6 * temp_original)
        
        # Bootstrap arrays
        cp_boot <- numeric(B)
        cpk_boot <- numeric(B)
        cpm_boot <- numeric(B)
        
        # Bootstrap loop
        for(b in 1:B) {
          boot_indices <- sample(1:n, n, replace = TRUE)
          sample <- sample_fixed[boot_indices]
          
          x_bar <- mean(sample)
          s <- sd(sample)
          temp <- sqrt(s^2 + (x_bar - T)^2)
          
          cp_boot[b] <- (USL - LSL)/(6 * s)
          cpk_boot[b] <- min((USL - x_bar)/(3 * s), (x_bar - LSL)/(3 * s))
          cpm_boot[b] <- (USL - LSL)/(6 * temp)
        }
        
        # Sort bootstrap estimates
        cp_boot <- sort(cp_boot)
        cpk_boot <- sort(cpk_boot)
        cpm_boot <- sort(cpm_boot)
        
        # Standard Bootstrap limits
        cp_boot_std <- sd(cp_boot)
        cpk_boot_std <- sd(cpk_boot)
        cpm_boot_std <- sd(cpm_boot)
        
        SB_L_cp <- cp_original - z_value * cp_boot_std
        SB_L_cpk <- cpk_original - z_value * cpk_boot_std
        SB_L_cpm <- cpm_original - z_value * cpm_boot_std
        
        # Percentile Bootstrap limits
        PB_L_cp <- quantile(cp_boot, alpha)
        PB_L_cpk <- quantile(cpk_boot, alpha)
        PB_L_cpm <- quantile(cpm_boot, alpha)
        
        # Bias-Corrected Percentile Bootstrap
        P0_cp <- sum(cp_boot <= cp_original)/B
        P0_cpk <- sum(cpk_boot <= cpk_original)/B
        P0_cpm <- sum(cpm_boot <= cpm_original)/B
        
        Z0_cp <- qnorm(P0_cp)
        Z0_cpk <- qnorm(P0_cpk)
        Z0_cpm <- qnorm(P0_cpm)
        
        PL_cp <- pnorm(2 * Z0_cp - z_value)
        PL_cpk <- pnorm(2 * Z0_cpk - z_value)
        PL_cpm <- pnorm(2 * Z0_cpm - z_value)
        
        idx_cp <- max(1, round(PL_cp * B))
        idx_cpk <- max(1, round(PL_cpk * B))
        idx_cpm <- max(1, round(PL_cpm * B))
        
        BCPB_L_cp <- cp_boot[idx_cp]
        BCPB_L_cpk <- cpk_boot[idx_cpk]
        BCPB_L_cpm <- cpm_boot[idx_cpm]
        
        # Update counters
        sb_counts <- sb_counts + c(SB_L_cp < true_cp,
                                   SB_L_cpk < true_cpk,
                                   SB_L_cpm < true_cpm)
        
        pb_counts <- pb_counts + c(PB_L_cp < true_cp,
                                   PB_L_cpk < true_cpk,
                                   PB_L_cpm < true_cpm)
        
        bcpb_counts <- bcpb_counts + c(BCPB_L_cp < true_cp,
                                       BCPB_L_cpk < true_cpk,
                                       BCPB_L_cpm < true_cpm)
      }
      
      # Calculate coverage probabilities
      cp_prob_SB <- sb_counts/N
      cp_prob_PB <- pb_counts/N
      cp_prob_BCPB <- bcpb_counts/N
      
      # Display results
      cat(sprintf("mu: %.1f, sigma: %.1f, n: %d\n", mu, sigma, n))
      cat(sprintf("True Values - Cp: %.2f, Cpk: %.2f, Cpm: %.2f\n",
                  true_cp, true_cpk, true_cpm))
      cat("Coverage Probabilities:\n")
      cat("        Cp      Cpk     Cpm\n")
      cat(sprintf("SB:     %.3f   %.3f   %.3f\n", cp_prob_SB[1], cp_prob_SB[2], cp_prob_SB[3]))
      cat(sprintf("PB:     %.3f   %.3f   %.3f\n", cp_prob_PB[1], cp_prob_PB[2], cp_prob_PB[3]))
      cat(sprintf("BCPB:   %.3f   %.3f   %.3f\n\n", cp_prob_BCPB[1], cp_prob_BCPB[2], cp_prob_BCPB[3]))
    }
  }
}