# Implementation of Various Process Capability Indices (PCIs)

A simulation study and implementation of various process capability indices (PCIs) confidence interval calculation methods.

## Overview

This project implements and compares different methods for calculating confidence intervals for process capability indices, including:
- Generalized confidence interval (Gpk)
- Bissell's approximation (Bpk)
- Kushler and Hurley's approximation (KHpk)
- Heavlin's approximation (Hpk)
- Nagata and Nagahata's approximation (Npk)

## Features

- Calculates and compares coverage probabilities for different methods
- Supports multiple process capability indices (Cpk, Cpmk, C''pk)
- Interactive visualization of results
- Parallel computation support for simulation studies
- Comprehensive statistical analysis tools

## Requirements

### R Dependencies
```R
library(parallel)
library(knitr)
library(kableExtra)
library(Rcpp)
```

### MATLAB Dependencies
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (optional, for improved performance)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/process-capability-analysis.git
cd process-capability-analysis
```

2. Install R dependencies:
```R
install.packages(c("parallel", "knitr", "kableExtra", "Rcpp"))
```

## Usage

### R Implementation

```R
# Run the main simulation
source("paper1.R")

# For custom parameters:
results <- simulate_Cpk(
    n = 30,                    # Sample size
    mu = 10,                   # Process mean
    sigma = 1,                 # Process standard deviation
    L = 7,                     # Lower specification limit
    U = 14,                    # Upper specification limit
    Cpk = 1.33,               # Target Cpk value
    alpha = 0.05,             # Significance level
    num_outer_simulations = 10000,
    num_inner_simulations = 10000
)
```

### MATLAB Implementation

```matlab
% Run the main simulation
simulateCpkTable()

% For custom analysis
results = simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, ...
    num_outer_simulations, num_inner_simulations);
```

## Key Functions

- `calculate_Cpk`: Calculates process capability index
- `calculate_Gpk_cpp`: C++ implementation of Gpk calculation (via Rcpp)
- `simulate_Cpk`: Performs Monte Carlo simulation
- `format_results`: Formats and displays results in tables
- `calculate_other_limits`: Computes various confidence interval approximations

## Results Interpretation

The simulation results include:
- Coverage probabilities for each method
- Expected values of confidence limits
- Comparative analysis across different sample sizes
- Performance evaluation at different confidence levels

## License

This project is available for academic and educational purposes.

## References

Mathew, T., Sebastian, G., & Kurian, K. M. (2007). Generalized confidence intervals for process capability indices. Quality and Reliability Engineering International, 23(4), 471-481.
