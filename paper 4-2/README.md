# Process Capability Indices Calculator
Implementation of algorithms for calculating lower confidence bounds of CPU and CPL process capability indices, as described in [Pearn et al. (2003)](https://doi.org/10.1016/S0026-2714(02)00264-0).

## Features
- Calculates lower confidence bounds for CPU and CPL indices
- Generates comprehensive parameter tables for different sample sizes
- Computes process capability indices using UMVUE estimators
- Supports multiple confidence levels
- Handles various sample sizes from n=5 to n=200
- Provides numerical solutions for practical implementation
- Includes PPM (parts per million) conversion tables

## Installation
```bash
git clone https://github.com/yourusername/process-capability-calculator
cd process-capability-calculator
```

## Usage
### MATLAB
```matlab
% Calculate lower confidence bounds for a dataset
lcb(80, 650, 0.95)

% Generate PPM conversion table
table1

% Generate lower confidence bounds table
table3

% Generate comparison table with Chou's method
table4
```

## File Structure
- `lcb.m`: Main implementation for calculating lower confidence bounds
- `table1.m`: Generates PPM conversion table for different CPU values
- `table3.m`: Generates comprehensive lower confidence bounds table
- `table4.m`: Implements comparison between new method and Chou's method
- `Appendix.m`: Contains sample dataset and implementation example

## Requirements
### MATLAB
- Statistics and Machine Learning Toolbox (for norminv, normcdf functions)
- MATLAB R2019b or newer recommended

## Mathematical Components
The implementation includes several key mathematical components:
- Noncentral t-distribution calculations
- Gamma function computations
- Polynomial root finding
- Iterative numerical methods for bound determination
- UMVUE (Uniformly Minimum-Variance Unbiased Estimator) calculations

## Implementation Notes
- Correction factor bn-1 uses logarithmic calculations for numerical stability
- Iterative tolerance is set to 0.01
- Maximum iterations are set to 500
- Confidence level typically set to 0.95 (95%)
- All calculations use double precision

## License
This project is available for academic and educational purposes.

## References
1. Pearn, W.L., & Shu, M.-H. (2003). An algorithm for calculating the lower confidence bounds of CPU and CPL with application to low-drop-out linear regulators. Microelectronics Reliability, 43(3), 495-502.

2. Pearn, W.L., & Shu, M.-H. (2003). Erratum to "An algorithm for calculating the lower confidence bounds of CPU and CPL with application to low-drop-out linear regulators". Microelectronics Reliability, 43(8), 1349.
