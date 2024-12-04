library(stats)

# 定義積分函數
# 參數說明：
# C: 待求解的參數
# n: 樣本數
# Cpk_hat: Cpk 的估計值
# xi_hat: xi 的估計值
# gamma: 置信水平
integrate_func <- function(C, n, Cpk_hat, xi_hat, gamma) {
  # 計算積分上限參數
  b <- 3*C + abs(xi_hat)
  
  # 定義被積函數
  integrand <- function(t) {
    # 計算卡方分布的累積分布函數值
    G_val <- pchisq((n-1)*(b*sqrt(n)-t)^2 / (9*n*Cpk_hat^2), df=n-1)
    # 計算正態分布的密度函數值
    phi_val <- dnorm(t + xi_hat*sqrt(n)) + dnorm(t - xi_hat*sqrt(n))
    return(G_val * phi_val)
  }
  
  # 進行數值積分
  result <- integrate(integrand, lower=0, upper=b*sqrt(n))
  # 返回與目標置信水平的差值
  return(result$value - (1 - gamma))
}

# 計算置信下界的函數
calculate_lcb <- function(n, Cpk_hat, xi_hat=1.0, gamma=0.95) {
  # 使用單調搜索找到積分函數為 0 的點
  result <- uniroot(integrate_func, c(0, Cpk_hat), 
                    n=n, Cpk_hat=Cpk_hat, xi_hat=xi_hat, gamma=gamma)
  return(result$root)
}

# 開始計時
start_time <- Sys.time()

# 設定參數範圍
n_values <- seq(10, 200, by = 5)  # 樣本數範圍：從 10 到 200，間隔 5
c_pk_hat_values_A <- seq(0.7, 1.8, by = 0.1)  # Cpk 估計值範圍 A ： 0.7 到 1.8
c_pk_hat_values_B <- seq(1.9, 3.0, by = 0.1)  # Cpk 估計值範圍 B ： 1.9 到 3.0

# 計算 Panel A 的結果
# 對每個 Cpk 值和樣本量組合計算置信下界
panel_A <- sapply(c_pk_hat_values_A, function(c_pk) {
  sapply(n_values, function(n) round(calculate_lcb(n, c_pk), 3))
})

# 計算 Panel B 的結果
panel_B <- sapply(c_pk_hat_values_B, function(c_pk) {
  sapply(n_values, function(n) round(calculate_lcb(n, c_pk), 3))
})

# 創建數據框
result_df_A <- data.frame(n = n_values)
result_df_A <- cbind(result_df_A, panel_A)
result_df_B <- data.frame(n = n_values)
result_df_B <- cbind(result_df_B, panel_B)

# 設定列名
colnames(result_df_A) <- c("n", as.character(c_pk_hat_values_A))
colnames(result_df_B) <- c("n", as.character(c_pk_hat_values_B))

# 輸出結果表格
cat("Panel A:\n")
print(result_df_A, row.names = FALSE)
cat("\nPanel B:\n")
print(result_df_B, row.names = FALSE)

# 計算並輸出總執行時間
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")
print(total_time)