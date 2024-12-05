# 定義函數計算給定 Cpk 值的 PPM 上下界
calculate_ppm_bounds <- function(cpk) {
  # 單邊規格的情況: 使用常態分配計算超出規格的 PPM 值
  upper_ppm <- (1 - pnorm(3 * cpk)) * 1e6
  
  # 雙邊規格的情況: 當製程置中時(對稱情況)的下界
  lower_ppm <- 2 * (1 - pnorm(3 * cpk)) * 1e6
  
  return(c(lower_ppm, upper_ppm))
}

# Cpk 從 0.6 到 2.0,間隔 0.01
cpk_seq <- seq(0.6, 2.0, by = 0.01)

# 對每個 Cpk 值計算其 PPM 上下界
results <- t(sapply(cpk_seq, calculate_ppm_bounds))

par(mar = c(5, 5, 2, 2))  # 調整邊界大小

# 建立主要圖形
plot(cpk_seq, results[,1], type = "l", 
     xlab = expression(C[pk]),    
     ylab = "PPM",              
     ylim = c(0, 80000),         
     xlim = c(0.6, 2.0),         
     yaxt = "n",                 
     cex.lab = 1.2,             # 調整標籤大小
     lwd = 1)                   # 設定線條寬度

# 加入上界線
lines(cpk_seq, results[,2], lwd = 1)

# 建立自訂的 Y 軸刻度
axis(2, at = seq(0, 80000, by = 20000), las = 1)

# 在圖形下方加入標題
title(sub = expression(paste("Fig. 1. The bounds on nonconforming units in PPM versus ", C[pk], ".")),
      line = 3.5)