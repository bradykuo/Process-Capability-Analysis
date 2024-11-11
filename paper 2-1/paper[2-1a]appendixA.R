library(stats)      # 用於統計計算
library(ggplot2)    # 用於繪圖
library(gridExtra)  # 用於排列多個圖表

# 定義計算p值的函數
calculate_p_value <- function(C, n_val, c_star, Cp) {
  # 定義G函數(卡方分配的累積分布函數)
  G <- function(c1, y) {
    pchisq((n_val - 1) * (3 * Cp * sqrt(n_val) - y)^2 / (9 * n_val * c1^2), df = n_val-1)
  }
  
  # 定義f函數(常態分配的機率密度函數)
  f <- function(y) {
    dnorm(y + 3 * (Cp - C) * sqrt(n_val)) + dnorm(y - 3 * (Cp - C) * sqrt(n_val))
  }
  
  # 定義被積函數
  integrand <- function(y) {
    G(c_star, y) * f(y)
  }
  
  # 設定積分上限
  upper_limit <- 3 * Cp * sqrt(n_val)
  # 執行數值積分
  result <- integrate(integrand, lower = 0, upper = upper_limit)
  
  return(result$value)
}

# 設定輸入參數
C <- 1.00  # 製程能力要求
n_values <- c(10, 30, 50, 100, 150, 200)  # 樣本數大小
c_star <- 1.15  # 臨界值
Cp_values <- seq(1, 1.9, by = 0.1)  # Cp值範圍

# 創建列表存儲所有圖表
plots_list <- list()

# 定義x軸刻度
x_breaks <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)

# 對每個樣本數計算p值並繪圖
for (n_val in n_values) {
  # 計算當前樣本數的p值
  p_values <- numeric(length(Cp_values))
  
  # 使用try-catch處理可能的錯誤
  for (i in seq_along(Cp_values)) {
    tryCatch({
      p_values[i] <- calculate_p_value(C, n_val, c_star, Cp_values[i])
    }, error = function(e) {
      cat("Error for n =", n_val, "and Cp =", Cp_values[i], ":", e$message, "\n")
      p_values[i] <- NA
    })
  }
  
  # 創建繪圖用的數據框
  plot_data <- data.frame(
    p_value = p_values,
    Cp = Cp_values
  )
  
  # 創建個別圖表
  p <- ggplot(plot_data, aes(x = Cp, y = p_value)) +
    geom_line(color = "blue", size = 1) +  # 繪製藍色線條
    theme_minimal() +  # 使用簡潔主題
    labs(
      title = paste("n =", n_val),
      x = "Cp",
      y = "P-value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    ) +
    scale_x_continuous(
      limits = c(1, 1.9),
      breaks = x_breaks
    )
  
  plots_list[[length(plots_list) + 1]] <- p
}

# 將所有圖表排列在一個網格中
grid.arrange(
  grobs = plots_list,
  ncol = 2,  
  nrow = 3,  
  top = "P-value vs Cp for Different Sample Sizes"
)