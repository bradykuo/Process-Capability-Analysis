# 定義製程參數
USL <- 61  # 規格上限
LSL <- 40  # 規格下限
T <- 49    # 目標值

# 定義不同的平均值和標準差組合
mu <- c(50, 52, 50, 52, 50, 52)      
sigma <- c(2.0, 2.0, 3.0, 3.0, 3.7, 3.7)  

# 初始化製程能力指數陣列
Cp <- numeric(length(mu))    
Cpk <- numeric(length(mu))   
Cpm <- numeric(length(mu))  

# 計算每種組合的製程能力指數
for(i in 1:length(mu)) {
  # 計算Cp (衡量製程的潛在能力)
  Cp[i] <- round((USL - LSL)/(6*sigma[i]), 2)
  
  # 計算Cpk (考慮製程的偏移程度)
  Cpu <- (USL - mu[i])/(3*sigma[i])  # 上限能力指數
  Cpl <- (mu[i] - LSL)/(3*sigma[i])  # 下限能力指數
  Cpk[i] <- round(min(Cpu, Cpl), 2)  # 取較小值
  
  # 計算Cpm (考慮製程的目標偏移)
  sigma_squared <- sigma[i]^2 + (mu[i] - T)^2  # 考慮變異與偏移
  Cpm[i] <- round((USL - LSL)/(6*sqrt(sigma_squared)), 2)
}

# 建立資料框
data <- data.frame(
  mu = mu,
  sigma = sigma,
  Cp = Cp,
  Cpk = Cpk,
  Cpm = Cpm
)

# 重新命名欄位
colnames(data) <- c("μ", "σ", "Cp", "Cpk", "Cpm")

# 顯示表格
print(data)