# Manufacturing Capability Control for Multiple Power-Distribution Switch Processes Based on Modified Cpk MPPAC

An implementation in R and MATLAB of manufacturing capability control methods for analyzing multiple power-distribution switch processes, based on the research paper by W.L. Pearn and Ming-Hung Shu (2003) published in Microelectronics Reliability.

## Overview

This package provides tools for:
- Calculating process capability indices (Cpk)
- Computing lower confidence bounds
- Generating Modified Process Performance Analysis Charts (MPPAC)
- Analyzing multiple process performance metrics

## Installation

### R Requirements
```R
install.packages(c("stats", "ggplot2", "gridExtra", "grid"))
```

### MATLAB Requirements
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

## Usage

### R Implementation

#### Calculating Lower Confidence Bounds
```R
# Calculate lower confidence bound
calculate_lcb <- function(n, Cpk_hat, xi_hat=1.0, gamma=0.95) {
    result <- uniroot(integrate_func, c(0, Cpk_hat), 
                     n=n, Cpk_hat=Cpk_hat, xi_hat=xi_hat, gamma=gamma)
    return(result$root)
}

# Example usage
result <- calculate_lcb(n=100, Cpk_hat=1.5)
```

#### Creating MPPAC Plots
```R
# Generate plots for different Cpk values
Cpk_hat_values <- c(0.7, 0.9, 1.2, 2.0, 2.5, 3.0)
n_values <- c(30, 50, 70, 100, 150, 200)
xi_values <- seq(0, 3, by = 0.1)

# Create plots using provided graph.R script
source("paper[2-2]graph.R")
```

### MATLAB Implementation

#### Calculating Lower Confidence Bounds
```matlab
% Calculate lower confidence bound
function lcb = calculate_lcb(n, Cpk_hat, xi_hat, gamma)
    if nargin < 3
        xi_hat = 1.0;
    end
    if nargin < 4
        gamma = 0.95;
    end
    
    options = optimset('Display', 'off');
    lcb = fzero(@(C) integrate_func(C, n, Cpk_hat, xi_hat, gamma), [0, Cpk_hat], options);
end

% Example usage
lcb = calculate_lcb(100, 1.5);
```

#### Creating MPPAC Plots
```matlab
% Run the MPPAC analysis
plotCpkAnalysis()
```

## Features

### Process Capability Analysis
- Calculation of Cpk index
- Lower confidence bound computation
- Integration with confidence levels
- Multiple process performance visualization

### MPPAC Visualization
- Capability zones representation
- Multiple process comparison
- Target line visualization
- Confidence bound plotting

### Process Categories
The analysis categorizes processes into:
- Super: Cpk ≥ 2.00 (NC < 0.002 ppm)
- Excellent: 1.67 ≤ Cpk < 2.00 (NC < 0.54 ppm)
- Satisfactory: 1.33 ≤ Cpk < 1.67 (NC < 66 ppm)
- Capable: 1.00 ≤ Cpk < 1.33 (NC < 2700 ppm)
- Inadequate: Cpk < 1.00 (NC > 2700 ppm)

Where NC represents non-conformities in parts per million (ppm).

## Mathematical Background

Key formulas implemented:

1. Process Capability Index (Cpk):
```
Cpk = min((USL - μ)/(3σ), (μ - LSL)/(3σ))
```

2. Lower Confidence Bound Integration:
```
∫[0 to b√n] G((n-1)(b√n-t)²/(9nx²)) * [φ(t + ξ√n) + φ(t - ξ√n)]dt
```
where:
- G() is the chi-square CDF
- φ() is the normal PDF
- n is sample size
- ξ is the standardized mean deviation

## File Structure

- `paper[2-2]graph.R`: R implementation of MPPAC plotting
- `paper[2-2]table.R`: R implementation of confidence bound tables
- `paper2_2graph.m`: MATLAB implementation of MPPAC plotting
- `paper2_2table.m`: MATLAB implementation of confidence bound tables

## License

This project is available for academic and educational purposes.

## References

1. Pearn, W.L., & Shu, M.-H. (2003). Manufacturing capability control for multiple power-distribution switch processes based on modified Cpk MPPAC. Microelectronics Reliability, 43(7), 963-975.
