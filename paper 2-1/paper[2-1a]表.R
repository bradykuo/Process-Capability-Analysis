# 定義計算臨界值c0的函數
calculate_c0 <- function(C, Cp, n, alpha) {
  # 定義被積分函數
  # G_term 是卡方分配的累積分配函數
  # f_term 是折疊常態分配的密度函數
  integrand <- function(y, c0) {
    G_term <- pchisq((n - 1) * (3 * Cp * sqrt(n) - y)^2 / (9 * n * c0^2), df = n - 1)
    f_term <- dnorm(y + 3 * (Cp - C) * sqrt(n)) + dnorm(y - 3 * (Cp - C) * sqrt(n))
    G_term * f_term
  }
  
  # 定義目標函數: 積分結果減去α值
  objective <- function(c0) {
    integrate(function(y) integrand(y, c0), lower = 0, upper = 3 * Cp * sqrt(n))$value - alpha
  }
  
  # 尋找適當的區間進行根搜尋
  # 初始值設定在0.5到5之間
  lower <- 0.5
  upper <- 5
  while (objective(lower) * objective(upper) > 0) {
    if (objective(lower) > 0) {
      lower <- lower / 2
    } else {
      upper <- upper * 2
    }
  }
  
  # 使用uniroot找出臨界值c0
  uniroot(objective, interval = c(lower, upper), tol = 1e-8)$root
}

# 記錄開始時間
start_time <- Sys.time()

# 設定參數值範圍
C_values <- c(1.00, 1.33, 1.5, 1.67, 2.00)  # Cpk要求值
n_values <- seq(10, 405, by = 5)  # 樣本數範圍，從10到405，每隔5取一值
alpha_values <- c(0.01, 0.025, 0.05)  # 顯著水準

# 對每個C值進行迴圈計算並生成表格
for (C in C_values) {
  # 建立矩陣儲存結果
  results <- matrix(nrow = length(n_values), ncol = length(alpha_values))
  
  # 對n值和α值進行雙重迴圈計算c0
  for (i in seq_along(n_values)) {
    n <- n_values[i]
    # 根據n的大小決定Cp值
    Cp <- if (n < 100) C + 0.33 else C + 0.12
    for (j in seq_along(alpha_values)) {
      alpha <- alpha_values[j]
      results[i, j] <- calculate_c0(C, Cp, n, alpha)
    }
  }
  
  # 將結果轉換為數據框
  table_df <- data.frame(n = n_values,
                         `α = 0.01` = results[, 1],
                         `α = 0.025` = results[, 2],
                         `α = 0.05` = results[, 3])
  
  # 使用kable套件顯示表格
  library(knitr)
  print(kable(table_df, col.names = c("n", "α = 0.01", "α = 0.025", "α = 0.05"),
              caption = paste("C =", C, "的臨界值c0表，n = 10(5)230，α = 0.01, 0.025, 0.05"),
              align = 'r', digits = 3))
}

# 計算並顯示總執行時間
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")
print(total_time)