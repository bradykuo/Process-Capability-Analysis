library(parallel)     # 用於平行運算
library(knitr)       # 用於產生報表
library(kableExtra)  # 用於美化表格
library(Rcpp)        # 用於整合 C++ 程式碼

# 使用 C++ 實作 Gpk 計算，以提升效能
cppFunction('
NumericVector calculate_Gpk_cpp(double x_bar, double s2, int n, double L, double U, int num_simulations) {
  double d = (U - L) / 2;           // 計算規格半寬
  double M = (U + L) / 2;           // 計算規格中點
  
  NumericVector Z = rnorm(num_simulations);        // 產生標準常態分配的隨機數
  NumericVector U2 = rchisq(num_simulations, n-1); // 產生卡方分配的隨機數
  
  // 計算 T 統計量
  NumericVector T_mu = x_bar - sqrt((n-1.0)/n) * (Z * sqrt(s2) / sqrt(n));
  NumericVector T_sigma2 = s2 * (n-1) / U2;
  
  // 計算並回傳 Gpk 值
  return (d - abs(T_mu - M)) / (3 * sqrt(T_sigma2));
}
')

# 計算 Cpk 值的函數
calculate_Cpk <- function(mu, sigma, L, U) {
  min((U - mu) / (3 * sigma), (mu - L) / (3 * sigma))
}

# 計算其他能力指標界限的函數
calculate_other_limits <- function(x_bar, s, n, L, U, alpha) {
  Cpk_hat <- calculate_Cpk(x_bar, s, L, U)
  z <- qnorm(1 - alpha)  # 計算標準常態分配的分位數
  
  # 將公式 key 入
  Bpk <- Cpk_hat - z * sqrt(1 / (9 * n) + Cpk_hat^2 / (2 * (n - 1)))
  KHpk <- Cpk_hat * (1 - z / sqrt(2 * (n - 1)))
  Hpk <- Cpk_hat - z * sqrt((n - 1) / (9 * n * (n - 3)) + Cpk_hat^2 * (1 / (2 * (n - 3)) * (1 + 6 / (n - 1))))
  Npk <- sqrt(1 - 2 / (5 * (n - 1))) * Cpk_hat - z * sqrt(Cpk_hat^2 / (2 * (n - 1)) + 1 / (9 * n))
  
  c(Bpk, KHpk, Hpk, Npk)
}

# 模擬計算 Cpk 的函數
simulate_Cpk <- function(n, mu, sigma, L, U, Cpk, alpha, num_outer_simulations, num_inner_simulations) {
  true_Cpk <- Cpk
  
  # 生成隨機樣本
  x <- matrix(rnorm(n * num_outer_simulations, mu, sigma), ncol = n)
  x_bar <- rowMeans(x)        # 計算樣本平均值
  s <- apply(x, 1, sd)        # 計算樣本標準差
  s2 <- s^2                   # 計算樣本變異數
  
  # 計算 Gpk 值
  Gpk_values <- sapply(seq_len(num_outer_simulations), function(i) {
    quantile(calculate_Gpk_cpp(x_bar[i], s2[i], n, L, U, num_inner_simulations), alpha)
  })
  
  # 計算其他指標值
  other_limits <- t(sapply(seq_len(num_outer_simulations), function(i) {
    calculate_other_limits(x_bar[i], s[i], n, L, U, alpha)
  }))
  
  limits <- cbind(Gpk_values, other_limits)
  
  # 計算覆蓋機率和期望值
  coverage <- colMeans(limits <= true_Cpk)
  expected_values <- colMeans(limits)
  
  list(coverage = round(coverage, 4), expected_values = round(expected_values, 4))
}

format_results <- function(results) {
  # 將巢狀列表轉換為資料框
  formatted_output <- do.call(rbind, lapply(names(results), function(Cpk) {
    do.call(rbind, lapply(names(results[[Cpk]]), function(n) {
      do.call(rbind, lapply(names(results[[Cpk]][[n]]), function(alpha) {
        c(Cpk, n, 1-as.numeric(alpha),
          results[[Cpk]][[n]][[alpha]]$coverage,
          results[[Cpk]][[n]][[alpha]]$expected_values)
      }))
    }))
  }))
  
  # 設定欄位名稱
  formatted_output <- as.data.frame(formatted_output)
  colnames(formatted_output) <- c("Cpk", "n", "1-α", 
                                  "Gpk", "Bpk", "KHpk", "Hpk", "Npk",
                                  "E(Gpk)", "E(Bpk)", "E(KHpk)", "E(Hpk)", "E(Npk)")
  
  # 轉換數值型態
  formatted_output$Cpk <- as.numeric(as.character(formatted_output$Cpk))
  formatted_output$n <- as.numeric(as.character(formatted_output$n))
  formatted_output$`1-α` <- as.numeric(as.character(formatted_output$`1-α`))
  
  # 四捨五入數值欄位
  numeric_cols <- 4:13
  formatted_output[numeric_cols] <- lapply(formatted_output[numeric_cols], 
                                           function(x) round(as.numeric(as.character(x)), 4))
  
  return(formatted_output)
}

# 設定隨機種子序，確保結果可重現
set.seed(123)  

# 設定模擬參數
Cpk_values <- c(1, 1.33, 1.5, 2, 2.5, 3)
n_values <- c(10, 20, 30, 40, 50)
alpha_values <- c(0.1, 0.05)
num_outer_simulations <- 10000
num_inner_simulations <- 10000
L <- 7    # 規格下限
U <- 14   # 規格上限
mu <- 10  # 製程平均值
d <- (U - L) / 2    # 規格半寬
M <- (U + L) / 2    # 規格中點

start_time <- Sys.time()  # 記錄開始時間

# 執行平行運算模擬
results <- mclapply(Cpk_values, function(Cpk) {
  sigma <- (d - abs(mu - M)) / (3 * Cpk)  # 計算對應的標準差
  
  lapply(n_values, function(n) {
    lapply(alpha_values, function(alpha) {
      cat(sprintf("Processing Cpk = %.2f, n = %d, alpha = %.2f\n", Cpk, n, alpha))
      simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, num_outer_simulations, num_inner_simulations)
    })
  })
}, mc.cores = detectCores())  # 使用所有可用的 CPU 核心

# 設定結果列表的名稱
names(results) <- as.character(Cpk_values)
for (i in seq_along(results)) {
  names(results[[i]]) <- as.character(n_values)
  for (j in seq_along(results[[i]])) {
    names(results[[i]][[j]]) <- as.character(alpha_values)
  }
}

end_time <- Sys.time()  # 記錄結束時間
total_time <- difftime(end_time, start_time, units = "secs")

# 格式化並顯示結果
formatted_results <- format_results(results)

# 修改欄位名稱以增加間距
names(formatted_results)[3] <- "1-α   "
names(formatted_results)[8] <- "Npk   "

# 使用 kable 產生美化的表格
print(kable(formatted_results, 
            format = "pipe",
            digits = 4,
            align = c('r', 'r', 'r', rep('r', 5), rep('r', 5))) %>%
        add_header_above(c(" " = 3, 
                           "Coverage probability" = 5,
                           "Expected value" = 5)) %>%
        add_header_above(c(" " = 13)) %>%
        kable_styling(full_width = F))

# 印出執行時間
cat(sprintf("\nTotal execution time: %.2f 秒\n", total_time))