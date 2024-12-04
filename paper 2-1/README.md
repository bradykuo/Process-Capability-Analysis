# Process Capability Analysis Tools

This project provides various numerical methods for calculating critical values in process capability analysis using the Cpk index. It includes implementations in both MATLAB and R for robust industrial process capability testing.

## Overview

The project implements different numerical approaches to solve for critical values (c0) in process capability analysis, offering multiple optimization and integration techniques to ensure reliable results across different parameter ranges.

## Available Implementations

### MATLAB
- `fminbnd_method.m`: Uses MATLAB's bounded minimization solver
- `fsolve_method.m`: Implements nonlinear equation solver
- `lsqnonlin_method.m`: Uses nonlinear least squares optimization
- `optimize_method.m`: Global optimization approach with GlobalSearch
- `paper2_1lineplot.m`: Line plot visualizations for critical value analysis
- `paper2_1a_surfaceplot.m`: 3D surface plot visualization for Cpk analysis across multiple parameters
- `paper2_1a_table.m`: Table generation for critical values

### R
- `multiroot.R`: Implementation using the rootSolve package
- `nleqslv.R`: Implementation using nonlinear equation solver
- `paper[2-1a]table.R`: Table generation for critical values
- `paper[2-1a]appendixA.R`: P-value calculations
- `paper[2-1a]appendixB.R`: Critical value calculations
- `paper[2-1a]line.R`: Line plot visualizations
- `paper[2-1a]surface.R`: Surface plot visualizations

## Features

- Critical value (c0) calculation using multiple numerical methods
- P-value computation for process capability assessment
- Visualization tools:
  - 3D surface plots for analyzing relationships between Cpk, n, and c0
  - Line plots for parameter relationships
  - Interactive visualizations in MATLAB and R
- Table generation for critical values
- Robust error handling and numerical stability measures

## Usage Examples

### MATLAB Surface Plot
```matlab
% Generate surface plots for different Cpk values
plot_cpk_analysis()  % Creates 5 subplots for Cpk values [1.00, 1.33, 1.50, 1.67, 2.00]
```

### MATLAB Critical Value Calculation
```matlab
% Calculate critical value using fminbnd
C = 1.33;
n = 50;
alpha = 0.05;
Cp = C + (n < 100)*0.33 + (n >= 100)*0.12;
c0 = calculate_c0_fminbnd(C, Cp, n, alpha);
```

### R Critical Value Calculation
```r
# Calculate critical value using nleqslv
C <- 1.33
n <- 50
alpha <- 0.05
Cp <- if (n < 100) C + 0.33 else C + 0.12
c0 <- calculate_c0_nleqslv(C, Cp, n, alpha)
```

## Requirements

### MATLAB
- MATLAB Optimization Toolbox
- MATLAB Statistics and Machine Learning Toolbox
- MATLAB Figure Graphics

### R
- rootSolve
- nleqslv
- stats
- ggplot2
- gridExtra

## Installation

1. Clone the repository
2. Install required dependencies for your preferred implementation
3. Run the example scripts to verify installation

## Visualization Capabilities

The project includes several visualization methods:
- 3D surface plots showing relationships between:
  - Sample size (n)
  - Process capability (Cp)
  - Critical values (c0)
- Line plots for different sample sizes
- Interactive plots with adjustable viewing angles
- Color-coded visualizations with colorbar indicators

## Implementation Details

### Critical Value Calculation Methods
- Bounded optimization (fminbnd)
- Nonlinear equation solving (fsolve)
- Least squares optimization (lsqnonlin)
- Global optimization
- Root finding methods (multiroot, nleqslv)

### Numerical Integration
Each implementation uses numerical integration techniques to evaluate the objective function, with proper error handling and boundary condition management.

## References

1. Testing process performance based on capability index Cpk with critical values (Computers & Industrial Engineering, 2004)
2. Process Capability Analysis Lecture 3-2: The Cpk CDF (Supplementary material)

## License

This project is available for academic and educational purposes.
