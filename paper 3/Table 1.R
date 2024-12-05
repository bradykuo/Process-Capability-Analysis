# 計算給定 Cpk 值的 PPM 界限
calculate_ppm_bounds <- function(cpk) {
  # 單邊規格: 計算超出規格上限的 PPM
  upper_ppm <- (1 - pnorm(3 * cpk)) * 1e6
  
  # 雙邊規格: 當製程居中時的下限 (對稱情況)
  lower_ppm <- 2 * (1 - pnorm(3 * cpk)) * 1e6
  
  return(c(lower_ppm, upper_ppm))
}

cpk_values <- c(0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.24, 1.25, 
                1.30, 1.33, 1.40, 1.45, 1.50, 1.60, 1.67, 1.70, 1.80, 1.90, 2.00)

# 計算每個 Cpk 值對應的 PPM 界限
results <- t(sapply(cpk_values, calculate_ppm_bounds))

cpk_table <- data.frame(
  Cpk = cpk_values,
  Lower_Bound = round(results[,1], 3),  # 下界 PPM
  Upper_Bound = round(results[,2], 3)   # 上界 PPM
)

# 輸出結果表格
print(cpk_table, row.names = FALSE)

# 驗證特定值與論文的對照
cat("\nVerification of selected values:\n")
cat("For Cpk = 1.33: Expected ≈ 66/33 PPM, Calculated:", 
    round(calculate_ppm_bounds(1.33), 3), "PPM\n")
cat("For Cpk = 1.67: Expected ≈ 0.544/0.272 PPM, Calculated:", 
    round(calculate_ppm_bounds(1.67), 3), "PPM\n")