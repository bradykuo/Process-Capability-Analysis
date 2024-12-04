# 計算臨界值的函數
calculate_critical_value <- function(C = 1.00, n = 38, alpha = 0.05, r_initial = 1.260) {
  # 根據樣本數設定 Cp
  # n < 100 時，Cp = C + 0.33
  # n >= 100 時，Cp = C + 0.12
  if (n < 100) {
    Cp <- C + 0.33
  } else {
    Cp <- C + 0.12
  }
  
  # 定義 G 函數 (卡方分配的累積分布函數)
  G <- function(c0, y) {
    pchisq((n - 1) * (3 * Cp * n^0.5 - y)^2 / (9 * n * c0^2), df = n-1)
  }
  
  # 定義 f 函數 (常態分配的機率密度函數)
  f <- function(y) {
    dnorm(y + 3 * (Cp - C) * n^0.5) + dnorm(y - 3 * (Cp - C) * n^0.5)
  }
  
  # 定義 pV 函數 (執行數值積分)
  pV <- function(c0) {
    # 計算積分上限
    upper_limit <- 3 * Cp * n^0.5
    # 定義被積函數
    integrand <- function(y) G(c0, y) * f(y)
    # 執行數值積分
    result <- integrate(integrand, lower = 0, upper = upper_limit)$value
    return(result)
  }
  
  # 臨界值計算函數 (遞迴搜尋)
  cv <- function(r) {
    c0 <- r
    p_val <- pV(c0)
    p_val_minus <- pV(c0 - 0.001)
    
    # 判斷是否找到臨界值
    if (p_val <= alpha && p_val_minus > alpha) {
      return(c0)
      # 如果 p 值太小，減小 c0
    } else if (p_val < alpha && p_val_minus < alpha) {
      return(cv(c0 - 0.001))
      # 如果 p 值太大，增大 c0
    } else {
      return(cv(c0 + 0.001))
    }
  }
  
  # 使用錯誤處理機制計算臨界值
  tryCatch({
    critical_value <- cv(r_initial)
    # 回傳臨界值和相關參數
    return(list(
      critical_value = critical_value,
      parameters = list(
        C = C,
        n = n,
        alpha = alpha,
        Cp = Cp
      )
    ))
  }, error = function(e) {
    return(paste("Error calculation:", e$message))
  })
}

# 使用範例值執行計算
result <- calculate_critical_value()
print(result$critical_value)

# 印出結果
cat("Parameter settings:\n")
cat(sprintf("Process capability requirements C = %.2f\n", result$parameters$C))
cat(sprintf("Sample number n = %d\n", result$parameters$n))
cat(sprintf("Significance level α = %.2f\n", result$parameters$alpha))
cat(sprintf("Process characteristic parameter Cp = %.2f\n", result$parameters$Cp))
cat(sprintf("\nCritical value = %.3f\n", result$critical_value))