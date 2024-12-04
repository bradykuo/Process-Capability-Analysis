library(stats)
library(ggplot2)
library(gridExtra)
library(grid)

# 定義積分函數，用於計算置信下界
# 參數說明：
# C: 待求解的參數
# n: 樣本數
# Cpk_hat: Cpk 的估計值
# xi_hat: xi 的估計值
# gamma: 置信水平
integrate_func <- function(C, n, Cpk_hat, xi_hat, gamma) {
  b <- 3*C + abs(xi_hat)
  
  # 定義被積函數
  # t: 積分變量
  # 返回：卡方分布和正態分布的組合
  integrand <- function(t) {
    G_val <- pchisq((n-1)*(b*sqrt(n)-t)^2 / (9*n*Cpk_hat^2), df=n-1)
    phi_val <- dnorm(t + xi_hat*sqrt(n)) + dnorm(t - xi_hat*sqrt(n))
    return(G_val * phi_val)
  }
  
  # 進行數值積分並返回與目標置信水平的差值
  result <- integrate(integrand, lower=0, upper=b*sqrt(n))
  return(result$value - (1 - gamma))
}

# 計算置信下界的函數
# 使用單調搜索(uniroot)找到積分函數為 0 的點
calculate_lcb <- function(n, Cpk_hat, xi_hat, gamma=0.95) {
  result <- uniroot(integrate_func, c(0, Cpk_hat), 
                    n=n, Cpk_hat=Cpk_hat, xi_hat=xi_hat, gamma=gamma)
  return(result$root)
}

# 計算不同 n 和 xi 值組合下的 C 值
calculate_C_vs_xi <- function(Cpk_hat, n_values, xi_values) {
  results <- expand.grid(n = n_values, xi = xi_values)
  results$C <- mapply(calculate_lcb, results$n, 
                      MoreArgs = list(Cpk_hat = Cpk_hat, gamma = 0.95), 
                      xi_hat = results$xi)
  return(results)
}

# 開始計時
start_time <- Sys.time()

# 設定參數範圍
n_values <- c(30, 50, 70, 100, 150, 200)  # 樣本數範圍
xi_values <- seq(0, 3, by = 0.1)          # xi 值範圍
Cpk_hat_values <- c(0.7, 0.9, 1.2, 2.0, 2.5, 3.0)  # Cpk 估計值範圍

# 創建空列表存放繪圖結果
plots <- list()

# 對每個 Cpk_hat 值進行循環計算和繪圖
for (i in 1:length(Cpk_hat_values)) {
  Cpk_hat <- Cpk_hat_values[i]
  results <- calculate_C_vs_xi(Cpk_hat, n_values, xi_values)
  
  # 創建單個圖形
  p <- ggplot(results, aes(x = xi, y = C, color = factor(n))) +
    geom_line() +
    labs(title = paste("Ĉpk =", Cpk_hat),
         x = "|ξ|", y = "C") +
    theme_minimal() +
    theme(legend.position = "none") +  # 移除個別圖例
    scale_color_discrete(name = "n")
  
  plots[[i]] <- p
}

# 創建共同圖例
legend_plot <- ggplot(results, aes(x = xi, y = C, color = factor(n))) +
  geom_line() +
  scale_color_discrete(name = "n") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 提取圖例
tmp <- ggplot_gtable(ggplot_build(legend_plot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

# 將所有圖形排列在網格中
arranged_plots <- gridExtra::grid.arrange(
  grobs = c(plots, list(legend)),
  layout_matrix = rbind(matrix(1:6, nrow = 2, byrow = TRUE),
                        c(NA, 7, NA)),
  top = grid::textGrob("C vs |ξ| for different Ĉpk values", 
                       gp = grid::gpar(fontsize = 16, font = 2))
)

# 儲存圖形
ggsave("Cpk_plots.png", arranged_plots, width = 15, height = 10, dpi = 300)

# 結束計時並輸出總執行時間
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "secs")
print(total_time)