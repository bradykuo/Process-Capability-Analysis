library(ggplot2)    
library(pracma)     # 數值計算
library(tidyr)      # 資料整理
library(dplyr)      # 資料處理
library(metR)       # 繪製等高線圖
library(gridExtra)  # 圖形排版

alpha <- 0.05    # 生產者風險
beta <- 0.1      # 消費者風險
C_AQL <- 1.33    # 可接受品質水準
C_LTPD <- 1.00   # 容忍品質水準
xi <- 1.0        # 根據論文設定 ξ = 1.0
b1 <- 3 * C_AQL + abs(xi)
b2 <- 3 * C_LTPD + abs(xi)

# 計算 S1 函數：生產者風險相關
calculate_S1 <- function(n, c0) {
  integrand <- function(t) {
    G_term <- pchisq((n - 1) * (b1 * sqrt(n) - t)^2 / (9 * n * c0^2), df = n-1)
    phi_term <- dnorm(t + xi*sqrt(n)) + dnorm(t - xi*sqrt(n))
    return(G_term * phi_term)
  }
  
  upper_limit <- b1 * sqrt(n)
  result <- try(integrate(integrand, lower = 0, upper = upper_limit)$value - (1 - alpha), silent = TRUE)
  
  if(inherits(result, "try-error")) {
    return(NA)
  }
  return(result)
}

# 計算 S2 函數：消費者風險相關
calculate_S2 <- function(n, c0) {
  integrand <- function(t) {
    G_term <- pchisq((n - 1) * (b2 * sqrt(n) - t)^2 / (9 * n * c0^2), df = n-1)
    phi_term <- dnorm(t + xi*sqrt(n)) + dnorm(t - xi*sqrt(n))
    return(G_term * phi_term)
  }
  
  upper_limit <- b2 * sqrt(n)
  result <- try(integrate(integrand, lower = 0, upper = upper_limit)$value - beta, silent = TRUE)
  
  if(inherits(result, "try-error")) {
    return(NA)
  }
  return(result)
}

# 尋找特定 n 值下使 S1 = 0 的 c0 值
find_c0_for_S1 <- function(n) {
  f <- function(c0) calculate_S1(n, c0)
  tryCatch({
    result <- uniroot(f, interval = c(0.8, 1.5), tol = 1e-6)
    return(result$root)
  }, error = function(e) NA)
}

# 尋找特定 n 值下使 S2 = 0 的 c0 值
find_c0_for_S2 <- function(n) {
  f <- function(c0) calculate_S2(n, c0)
  tryCatch({
    result <- uniroot(f, interval = c(0.8, 1.5), tol = 1e-6)
    return(result$root)
  }, error = function(e) NA)
}

# 尋找 S1 和 S2 曲線的交點
find_intersection_precise <- function() {
  # 最小化 S1=0 和 S2=0 時的 c0 值差異
  f <- function(n) {
    n <- ceiling(n)  # 確保 n 為整數
    c0_s1 <- find_c0_for_S1(n)
    c0_s2 <- find_c0_for_S2(n)
    if(is.na(c0_s1) || is.na(c0_s2)) return(NA)
    return(c0_s1 - c0_s2)
  }
  
  # 使用 uniroot 尋找交點
  result <- tryCatch({
    n_intersection <- uniroot(f, interval = c(10, 150), tol = 1e-6)
    n <- ceiling(n_intersection$root)
    c0 <- find_c0_for_S1(n)
    return(c(n = n, c0 = c0))
  }, error = function(e) NULL)
  
  return(result)
}

# 尋找精確交點
intersection <- find_intersection_precise()

n_values <- seq(10, 150, length.out = 75)   # 增加解析度
c0_values <- seq(0.8, 1.5, length.out = 75) 
grid <- expand.grid(n = n_values, c0 = c0_values)

# 計算每個網格點的 S1 和 S2 值
grid$S1 <- mapply(calculate_S1, grid$n, grid$c0)
grid$S2 <- mapply(calculate_S2, grid$n, grid$c0)

# 繪製 S1 等高線圖
contour_plot_S1 <- ggplot(grid, aes(x = n, y = c0, z = S1)) +
  geom_contour(aes(z = S1), 
               bins = 10,
               color = "blue") +
  geom_contour(aes(z = S1), 
               breaks = 0,
               color = "blue",
               size = 1) +
  geom_text_contour(aes(z = S1), 
                    skip = 2,
                    rotate = TRUE,
                    size = 3,
                    stroke = 0.2) +
  scale_x_continuous(limits = c(10, 150), 
                     breaks = seq(0, 150, by = 25)) +
  scale_y_continuous(limits = c(0.8, 1.5), 
                     breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "Contour Plot of S1(n, c0)",
    x = "n",
    y = "c0"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# 繪製 S2 等高線圖
contour_plot_S2 <- ggplot(grid, aes(x = n, y = c0, z = S2)) +
  geom_contour(aes(z = S2), 
               bins = 10,
               color = "red") +
  geom_contour(aes(z = S2), 
               breaks = 0,
               color = "red",
               size = 1) +
  geom_text_contour(aes(z = S2), 
                    skip = 2,
                    rotate = TRUE,
                    size = 3,
                    stroke = 0.2) +
  scale_x_continuous(limits = c(10, 150), 
                     breaks = seq(0, 150, by = 25)) +
  scale_y_continuous(limits = c(0.8, 1.5), 
                     breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "Contour Plot of S2(n, c0)",
    x = "n",
    y = "c0"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# 尋找交點函數
find_intersection <- function(grid) {
  zero_threshold <- 0.01
  intersection_points <- grid[abs(grid$S1) < zero_threshold & abs(grid$S2) < zero_threshold, ]
  
  if(nrow(intersection_points) > 0) {
    closest_point <- intersection_points[which.min(intersection_points$S1^2 + intersection_points$S2^2), ]
    closest_point$n <- ceiling(closest_point$n)
    c0_value <- find_c0_for_S1(closest_point$n)
    return(c(n = closest_point$n, c0 = c0_value))
  } else {
    return(NULL)
  }
}

# 尋找交點
intersection <- find_intersection(grid)

# 繪製合併等高線圖
combined_plot <- ggplot(grid) +
  # 加入藍色 S1 等高線
  geom_contour(aes(x = n, y = c0, z = S1), 
               color = "blue",
               linetype = "solid") +
  geom_text_contour(aes(x = n, y = c0, z = S1),
                    skip = 1,
                    rotate = TRUE,
                    size = 3,
                    color = "blue",
                    stroke = 0.2) +
  
  # 加入紅色 S2 等高線
  geom_contour(aes(x = n, y = c0, z = S2),
               color = "red",
               linetype = "dashed") +
  geom_text_contour(aes(x = n, y = c0, z = S2),
                    skip = 1,
                    rotate = TRUE,
                    size = 3,
                    color = "red",
                    stroke = 0.2) +
  
  # 加入零等高線以突顯交點
  geom_contour(aes(x = n, y = c0, z = S1),
               breaks = 0,
               color = "blue",
               size = 1) +
  geom_contour(aes(x = n, y = c0, z = S2),
               breaks = 0,
               color = "red",
               size = 1) +
  
  # 加入交點和標籤
  geom_point(data = data.frame(n = intersection["n"], c0 = intersection["c0"]),
             aes(x = n, y = c0),
             color = "black",
             size = 3) +
  annotate("text", 
           x = intersection["n"], 
           y = intersection["c0"] + 0.05,
           label = sprintf("(n, c0) = (%.2f, %.4f)", 
                           intersection["n"], 
                           intersection["c0"]),
           size = 3.5,
           hjust = 0) +
  
  # 自訂座標軸和標籤
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 50)) +
  scale_y_continuous(limits = c(0.8, 1.5), breaks = seq(0.8, 1.5, by = 0.1)) +
  labs(
    title = "Combined Contour Plot of S1(n, c0) and S2(n, c0)",
    subtitle = "S1: Blue (solid), S2: Red (dashed)",
    x = "n",
    y = "C0"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

# 輸出交點值
print(sprintf("Intersection point: n = %d, c0 = %.4f", 
              round(intersection[1]), intersection[2]))

# 將所有圖形排列在網格中
grid.arrange(
  contour_plot_S1, 
  contour_plot_S2, 
  combined_plot,
  layout_matrix = rbind(c(1,2), c(3,3)),
  heights = c(1, 1)
)