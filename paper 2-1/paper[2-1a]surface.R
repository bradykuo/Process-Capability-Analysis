library(plotly)  # 載入 3D 互動式繪圖套件

# 計算臨界值 c0 的函數
calculate_c0 <- function(C, Cp, n, alpha) {
  # 定義被積分函數
  integrand <- function(y, c0) {
    # 計算卡方分配項
    G_term <- pchisq((n - 1) * (3 * Cp * sqrt(n) - y)^2 / (9 * n * c0^2), df = n - 1)
    # 計算常態分配項
    f_term <- dnorm(y + 3 * (Cp - C) * sqrt(n)) + dnorm(y - 3 * (Cp - C) * sqrt(n))
    G_term * f_term
  }
  
  # 定義目標函數
  objective <- function(c0) {
    integrate(function(y) integrand(y, c0), lower = 0, upper = 3 * Cp * sqrt(n))$value - alpha
  }
  
  # 尋找根值的區間
  lower <- 0.5
  upper <- 5
  while (objective(lower) * objective(upper) > 0) {
    if (objective(lower) > 0) {
      lower <- lower / 2
    } else {
      upper <- upper * 2
    }
  }
  
  # 使用 uniroot 找出精確值
  uniroot(objective, interval = c(lower, upper), tol = 1e-8)$root
}

# 創建曲面圖數據的函數
create_surface_data <- function(C, Cp_min, Cp_max) {
  alpha <- 0.05
  
  Cp_values <- seq(Cp_min, Cp_max, length.out = 50)  # Cp 值序列
  n_values <- seq(30, 300, length.out = 50)          # 樣本數 n 序列
  
  # 計算每個網格點的 c0 值
  c0_values <- outer(Cp_values, n_values, 
                     Vectorize(function(Cp, n) calculate_c0(C, Cp, n, alpha)))
  
  # 回傳曲面圖所需數據
  list(
    x = n_values,
    y = Cp_values,
    z = c0_values,
    type = "surface",
    name = sprintf("Cpk = %.2f", C),
    scene = sprintf("scene%d", which(C_values == C))
  )
}

# 設定所有圖形的參數
C_values <- c(1.00, 1.33, 1.50, 1.67, 2.00)  # Cpk 值
Cp_ranges <- list(  # 每個圖的 Cp 範圍
  c(1.00, 2.00),
  c(1.33, 2.33),
  c(1.50, 2.50),
  c(1.67, 2.67),
  c(2.00, 3.00)
)

# 創建子圖
scenes <- list()
for(i in 1:5) {
  scene_name <- sprintf("scene%d", i)
  scenes[[scene_name]] <- list(
    domain = list(  # 設定子圖位置
      x = c((i-1)%%3/3, ((i-1)%%3+1)/3),
      y = c(floor((i-1)/3)/2, (floor((i-1)/3)+1)/2)
    ),
    xaxis = list(title = "n"),
    yaxis = list(title = "Cp", autorange = "reversed"),
    zaxis = list(title = "c0"),
    camera = list(  # 設定視角
      eye = list(x = 1.5, y = 1.5, z = 1.5)
    ),
    aspectmode = 'cube'  # 設定軸比例
  )
}

# 生成所有圖形的數據
plot_data <- lapply(seq_along(C_values), function(i) {
  create_surface_data(C_values[i], Cp_ranges[[i]][1], Cp_ranges[[i]][2])
})

# 創建組合圖
fig <- plot_ly() %>%
  layout(
    title = "Surface plots of c0 for different Cpk values (α = 0.05)",
    showlegend = TRUE,
    scene = scenes$scene1,
    scene2 = scenes$scene2,
    scene3 = scenes$scene3,
    scene4 = scenes$scene4,
    scene5 = scenes$scene5
  )

# 添加每個曲面到圖形中
for(i in seq_along(plot_data)) {
  fig <- add_trace(fig, 
                   x = plot_data[[i]]$x,
                   y = plot_data[[i]]$y,
                   z = plot_data[[i]]$z,
                   type = "surface",
                   name = sprintf("Cpk = %.2f", C_values[i]),
                   scene = sprintf("scene%d", i))
}

# 印出組合圖
fig