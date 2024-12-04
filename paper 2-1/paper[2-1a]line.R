library(ggplot2)  # 用於繪製圖形
library(dplyr)    # 用於數據處理
library(gridExtra)  # 用於排列多個圖形
library(grid)     # 用於添加標題等圖形元素

# 定義計算臨界值 c0 的函數
calculate_c0 <- function(C, Cp, n, alpha) {
  # 定義被積分函數
  # G_term是卡方分配的累積分配函數
  # f_term是折疊常態分配的密度函數
  integrand <- function(y, c0) {
    G_term <- pchisq((n - 1) * (3 * Cp * sqrt(n) - y)^2 / (9 * n * c0^2), df = n - 1)
    f_term <- dnorm(y + 3 * (Cp - C) * sqrt(n)) + dnorm(y - 3 * (Cp - C) * sqrt(n))
    G_term * f_term
  }
  
  # 定義目標函數：積分結果減去 α 值
  objective <- function(c0) {
    tryCatch({
      result <- integrate(function(y) integrand(y, c0), 
                          lower = 0, 
                          upper = 3 * Cp * sqrt(n))$value - alpha
      return(result)
    }, error = function(e) {
      return(NA)
    })
  }
  
  # 設定初始搜索區間
  lower <- 0.5
  upper <- 5
  
  # 設定最大嘗試次數
  max_attempts <- 10
  attempts <- 0
  
  # 調整搜索區間直到找到合適的範圍
  while (attempts < max_attempts) {
    obj_lower <- objective(lower)
    obj_upper <- objective(upper)
    
    if (!is.na(obj_lower) && !is.na(obj_upper) && obj_lower * obj_upper < 0) {
      break
    }
    
    # 根據函數值調整區間
    if (is.na(obj_lower) || (!is.na(obj_lower) && obj_lower > 0)) {
      lower <- lower / 2
    }
    if (is.na(obj_upper) || (!is.na(obj_upper) && obj_upper < 0)) {
      upper <- upper * 2
    }
    
    attempts <- attempts + 1
  }
  
  # 使用 uniroot 找出臨界值 c0
  result <- try(uniroot(objective, interval = c(lower, upper), tol = 1e-6)$root, silent = TRUE)
  if (inherits(result, "try-error")) {
    return(NA)
  }
  return(result)
}

# 為每個 Cpk 值創建圖形的函數
create_cpk_plot <- function(Cpk, alpha = 0.05) {
  # 設定樣本數 n 範圍
  n_values <- c(30, 50, 70, 100, 150, 200, 250, 300)
  
  # 根據 Cpk 創建 Cp 範圍
  Cp_values <- seq(Cpk, Cpk + 1, length.out = 20)
  
  # 創建所有組合
  plot_data <- expand.grid(n = n_values, Cp = Cp_values)
  
  # 計算每個組合的 c0 值
  plot_data$c0 <- mapply(function(n, Cp) calculate_c0(Cpk, Cp, n, alpha),
                         plot_data$n, plot_data$Cp)
  
  # 計算 y 軸範圍
  y_max <- max(plot_data$c0, na.rm = TRUE) + 0.01
  y_min <- min(plot_data$c0, na.rm = TRUE)
  
  # 創建圖形
  p <- ggplot(plot_data, aes(x = Cp, y = c0, color = factor(n))) +
    geom_line(linewidth = 0.5) +               # 繪製折線
    scale_color_discrete(name = "Sample size (n)") +  # 設定圖例
    labs(title = sprintf("Cpk = %.2f", Cpk),   # 設定標題和軸
         x = "Cp",
         y = "c0") +
    theme_minimal() +          
    theme(
      legend.position = "right",              # 圖例位置
      plot.title = element_text(hjust = 0.5, size = 12),  # 標題置中和大小
      axis.title = element_text(size = 10),   # 軸標籤大小
      legend.title = element_text(size = 10),  # 圖例標題大小
      legend.text = element_text(size = 8)     # 圖例文字大小
    ) +
    coord_cartesian(ylim = c(y_min, y_max), expand = FALSE)  # 設定 y 軸範圍
  
  return(p)
}

# 設定要繪製的 Cpk 值
Cpk_values <- c(1.00, 1.33, 1.50, 1.67, 2.00)
plots <- list()

# 創建所有圖形
for(i in seq_along(Cpk_values)) {
  plots[[i]] <- create_cpk_plot(Cpk_values[i])
}

# 添加總標題並排列圖形
title <- textGrob("Plots of c0 versus Cp for different Cpk values (α = 0.05)", 
                  gp = gpar(fontsize = 14))
# 印出個別圖形
grid.arrange(
  plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
  ncol = 2,
  top = title
)
