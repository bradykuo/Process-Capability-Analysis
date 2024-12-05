library(ggplot2)  # 繪圖
library(dplyr)    # 資料處理
library(pracma)   # 數值計算

# 計算方程式(9)和(10)中的積分
calculate_integral <- function(n, c0, C_value, xi) {
  b <- 3 * C_value + abs(xi)
  
  integrand <- function(t) {
    # G 是自由度為 n-1 的卡方分配累積函數
    G_term <- pchisq((n - 1) * (b * sqrt(n) - t)^2 / (9 * n * c0^2), df = n-1)
    # phi 是標準常態分配的機率密度函數
    phi_term <- dnorm(t + xi*sqrt(n)) + dnorm(t - xi*sqrt(n))
    return(G_term * phi_term)
  }
  
  result <- integrate(integrand, lower = 0, upper = b*sqrt(n))$value
  return(result)
}

# 尋找給定參數下的 c0 值
find_c0 <- function(n, C_AQL, C_LTPD, alpha, beta, xi) {
  objective <- function(c0) {
    eq1 <- calculate_integral(n, c0, C_AQL, xi) - (1 - alpha)
    eq2 <- calculate_integral(n, c0, C_LTPD, xi) - beta
    return(eq1^2 + eq2^2)
  }
  
  result <- optimize(objective, interval = c(0.8, 2.0))
  return(result$minimum)
}

xi_values <- seq(0, 2, by = 0.1)          # ξ 值
CAQL_values <- c(1.33, 1.50, 1.67, 2.00)  # 可接受品質水準值
CLTPD <- 1.00                             # 可容忍品質水準
alpha <- 0.05                             # 生產者風險
beta <- 0.05                              # 消費者風險

# 初始化資料框
sample_size_data <- data.frame()      # 樣本數
critical_value_data <- data.frame()    # 臨界值

# 計算每個 CAQL 和 ξ 組合的值
for(CAQL in CAQL_values) {
  for(xi in xi_values) {
    # 尋找滿足兩個方程式的最小樣本數
    n_min <- 10  # 起始樣本數
    while(TRUE) {
      c0 <- find_c0(n_min, CAQL, CLTPD, alpha, beta, xi)
      
      # 檢查是否滿足方程式
      eq1 <- calculate_integral(n_min, c0, CAQL, xi) >= (1 - alpha)
      eq2 <- calculate_integral(n_min, c0, CLTPD, xi) <= beta
      
      if(eq1 && eq2) break
      n_min <- n_min + 1
      if(n_min > 100) break  # 安全停止條件
    }
    
    # 儲存結果
    sample_size_data <- rbind(sample_size_data, 
                              data.frame(xi = xi, n = n_min, CAQL = as.factor(CAQL)))
    critical_value_data <- rbind(critical_value_data, 
                                 data.frame(xi = xi, c0 = c0, CAQL = as.factor(CAQL)))
  }
}

# 繪製圖a：樣本數 vs ξ
plot_a <- ggplot(sample_size_data, aes(x = xi, y = n, color = CAQL)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Required Sample Size vs ξ",
    x = "ξ",
    y = "Required sample size n"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80"),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# 繪製圖b：臨界值 vs ξ
plot_b <- ggplot(critical_value_data, aes(x = xi, y = c0, color = CAQL)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(1.0, 1.6)) +
  labs(
    title = "Critical Acceptance Value vs ξ",
    x = "ξ",
    y = "Critical value Co"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80"),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# 顯示圖形
print(plot_a)
print(plot_b)