# Variables sampling inspection scheme for resubmitted lots based on the process capability index Cpk
Implementation of a sampling inspection scheme for resubmitted lots based on the process capability index Cpk, as described in [Wu et al. (2012)](https://doi.org/10.1016/j.ejor.2011.09.042).

## Features
- Calculates optimal sample sizes and critical values for various risk levels
- Generates Operating Characteristic (OC) curves for different resubmission counts
- Computes Average Sample Number (ASN) curves
- Supports multiple CAQL/CLTPD combinations
- Handles resubmitted lots with m=1 to m=10 resubmissions
- Provides comprehensive parameter tables for practical implementation

## Installation
```bash
git clone https://github.com/yourusername/sampling-inspection-scheme
cd sampling-inspection-scheme
```

## Usage
### MATLAB
```matlab
% Generate Table 1 (m=2 resubmissions)
table1

% Generate Table 2 (m=3 resubmissions)
table2

% Generate Table 3 (varying m from 1-10)
table3

% Generate OC curves
fig1

% Generate ASN curves
fig2
```

## File Structure
- `table1.m`: Generates sampling plan parameters for m=2 resubmissions
- `table2.m`: Generates sampling plan parameters for m=3 resubmissions
- `table3.m`: Generates sampling plan parameters for m=1 to m=10
- `fig1.m`: Creates OC curves comparing different resubmission counts
- `fig2.m`: Creates ASN curves for different resubmission scenarios
- `analysis.m`: Basic implementation for parameter calculation

## Requirements
### MATLAB
- Statistics and Machine Learning Toolbox (for chi2cdf, normpdf functions)
- Optimization Toolbox (for fsolve function)
- MATLAB R2019b or newer recommended

## Mathematical Components
The implementation includes several key mathematical components:
- Numerical integration for acceptance probability calculation
- Non-linear equation solving for parameter optimization
- Chi-square and normal distribution computations
- Boundary coefficient calculations

## Implementation Notes
- Process shift parameter (ξ) is set to 1.0
- Function and step tolerances are set to 1e-8
- Maximum iterations for optimization are set to 1000
- Sample sizes are rounded up to ensure integer values

## License
This project is available for academic and educational purposes.

## References
Wu, C.-W., Aslam, M., & Jun, C.-H. (2012). Variables sampling inspection scheme for resubmitted lots based on the process capability index Cpk. European Journal of Operational Research, 217(3), 560-566.
