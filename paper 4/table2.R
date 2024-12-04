# 定義常數
USL <- 61        # 規格上限
LSL <- 40        # 規格下限
T <- 49          # 目標值
alpha <- 0.05    # 顯著水準

# 設定不同情境的參數值
mu_list <- c(50, 52)           # 平均值清單
sigma_list <- c(2.0, 3.0, 3.7) # 標準差清單
n_list <- c(20, 40, 70)        # 樣本數清單

B <- 1000    # Bootstrap 重抽樣次數
N <- 1000    # 模擬重複次數
z_value <- qnorm(1 - alpha)    # 計算常態分配的臨界值

start_time <- Sys.time()   # 記錄開始時間

# 主要模擬迴圈
for(mu_idx in 1:length(mu_list)) {
  mu <- mu_list[mu_idx]
  
  for(sigma_idx in 1:length(sigma_list)) {
    sigma <- sigma_list[sigma_idx]
    
    # 計算真實的製程能力指數
    true_cp <- (USL - LSL)/(6 * sigma)                                    # 計算 Cp 值 
    true_cpk <- min((USL - mu)/(3 * sigma), (mu - LSL)/(3 * sigma))      # 計算 Cpk 值 
    true_cpm <- (USL - LSL)/(6 * sqrt(sigma^2 + (mu - T)^2))            # 計算 Cpm 值 
    
    for(n_idx in 1:length(n_list)) {
      n <- n_list[n_idx]
      
      # 初始化計數器（用於記錄下限小於真實值的次數）
      sb_counts <- numeric(3)   # Standard Bootstrap 
      pb_counts <- numeric(3)   # Percentile Bootstrap 
      bcpb_counts <- numeric(3) # Biased Corrected Percentile Bootstrap
      
      for(sim in 1:N) {
        # 從常態分配生成原始樣本
        sample_fixed <- rnorm(n, mean = mu, sd = sigma)
        
        # 計算原始樣本的統計量
        x_bar_original <- mean(sample_fixed)
        s_original <- sd(sample_fixed)
        temp_original <- sqrt(s_original^2 + (x_bar_original - T)^2)
        
        # 計算原始樣本的製程能力指數
        cp_original <- (USL - LSL)/(6 * s_original)
        cpk_original <- min((USL - x_bar_original)/(3 * s_original),
                            (x_bar_original - LSL)/(3 * s_original))
        cpm_original <- (USL - LSL)/(6 * temp_original)
        
        # 初始化 Bootstrap 陣列
        cp_boot <- numeric(B)
        cpk_boot <- numeric(B)
        cpm_boot <- numeric(B)
        
        # Bootstrap 重抽樣迴圈
        for(b in 1:B) {
          boot_indices <- sample(1:n, n, replace = TRUE)
          sample <- sample_fixed[boot_indices]
          
          x_bar <- mean(sample)
          s <- sd(sample)
          temp <- sqrt(s^2 + (x_bar - T)^2)
          
          # 計算重抽樣後的製程能力指數
          cp_boot[b] <- (USL - LSL)/(6 * s)
          cpk_boot[b] <- min((USL - x_bar)/(3 * s), (x_bar - LSL)/(3 * s))
          cpm_boot[b] <- (USL - LSL)/(6 * temp)
        }
        
        # 排序 Bootstrap 估計值
        cp_boot <- sort(cp_boot)
        cpk_boot <- sort(cpk_boot)
        cpm_boot <- sort(cpm_boot)
        
        # Standard Bootstrap 信賴區間下限
        cp_boot_std <- sd(cp_boot)
        cpk_boot_std <- sd(cpk_boot)
        cpm_boot_std <- sd(cpm_boot)
        
        SB_L_cp <- cp_original - z_value * cp_boot_std
        SB_L_cpk <- cpk_original - z_value * cpk_boot_std
        SB_L_cpm <- cpm_original - z_value * cpm_boot_std
        
        # Percentile Bootstrap 信賴區間下限
        PB_L_cp <- quantile(cp_boot, alpha)
        PB_L_cpk <- quantile(cpk_boot, alpha)
        PB_L_cpm <- quantile(cpm_boot, alpha)
        
        # Biased Corrected Percentile Bootstrap
        # 計算原始估計值在 Bootstrap 分配中的位置
        P0_cp <- sum(cp_boot <= cp_original)/B
        P0_cpk <- sum(cpk_boot <= cpk_original)/B
        P0_cpm <- sum(cpm_boot <= cpm_original)/B
        
        # 轉換為 z 分數
        Z0_cp <- qnorm(P0_cp)
        Z0_cpk <- qnorm(P0_cpk)
        Z0_cpm <- qnorm(P0_cpm)
        
        # 計算偏差修正後的百分位
        PL_cp <- pnorm(2 * Z0_cp - z_value)
        PL_cpk <- pnorm(2 * Z0_cpk - z_value)
        PL_cpm <- pnorm(2 * Z0_cpm - z_value)
        
        # 找出對應的 Bootstrap 樣本位置
        idx_cp <- max(1, round(PL_cp * B))
        idx_cpk <- max(1, round(PL_cpk * B))
        idx_cpm <- max(1, round(PL_cpm * B))
        
        # 取得偏差修正後的信賴區間下限
        BCPB_L_cp <- cp_boot[idx_cp]
        BCPB_L_cpk <- cpk_boot[idx_cpk]
        BCPB_L_cpm <- cpm_boot[idx_cpm]
        
        # 更新計數：檢查各方法的下限是否小於真實值
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
      
      # 計算覆蓋機率
      cp_prob_SB <- sb_counts/N
      cp_prob_PB <- pb_counts/N
      cp_prob_BCPB <- bcpb_counts/N
      
      # 顯示結果
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

end_time <- Sys.time()         # 記錄結束時間
total_time <- difftime(end_time, start_time, units = "secs")  # 計算總執行時間

# 印出執行時間
cat(sprintf("\nTotal execution time: %.2f seconds\n", total_time))