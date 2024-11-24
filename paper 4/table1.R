# Define process parameters
USL <- 61  # Upper Specification Limit
LSL <- 40  # Lower Specification Limit
T <- 49    # Target value

# Define different combinations of mean and standard deviation
mu <- c(50, 52, 50, 52, 50, 52)
sigma <- c(2.0, 2.0, 3.0, 3.0, 3.7, 3.7)

# Initialize arrays for indices
Cp <- numeric(length(mu))
Cpk <- numeric(length(mu))
Cpm <- numeric(length(mu))

# Calculate indices for each combination
for(i in 1:length(mu)) {
  # Calculate Cp
  Cp[i] <- round((USL - LSL)/(6*sigma[i]), 2)
  
  # Calculate Cpk
  Cpu <- (USL - mu[i])/(3*sigma[i])
  Cpl <- (mu[i] - LSL)/(3*sigma[i])
  Cpk[i] <- round(min(Cpu, Cpl), 2)
  
  # Calculate Cpm
  sigma_squared <- sigma[i]^2 + (mu[i] - T)^2
  Cpm[i] <- round((USL - LSL)/(6*sqrt(sigma_squared)), 2)
}

# Create the data frame
data <- data.frame(
  mu = mu,
  sigma = sigma,
  Cp = Cp,
  Cpk = Cpk,
  Cpm = Cpm
)

# Rename columns
colnames(data) <- c("μ", "σ", "Cp", "Cpk", "Cpm")

# Display the table
print(data)