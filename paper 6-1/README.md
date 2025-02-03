# Repetitive Group Sampling Plan Implementation

Implementation of a repetitive group sampling (RGS) plan by variables inspection for product acceptance determination based on the process capability index Cpk, as described in [Wu et al. (2015)](https://doi.org/10.1504/EJIE.2015.069341).

## Features
- Calculates Operating Characteristic (OC) curves for single sampling and RGS plans
- Computes optimal sample sizes and critical values for various risk levels (α, β)
- Generates comprehensive parameter tables for different quality levels (CAQL, CLTPD)
- Calculates PPM bounds for different Cpk values
- Provides Average Sample Number (ASN) comparisons between sampling plans
- Supports multiple quality level combinations (1.33/1.00, 1.50/1.33, 1.67/1.33, 2.00/1.67)

## Installation
```bash
git clone https://github.com/yourusername/variables-rgs-inspection
cd variables-rgs-inspection
```

## Usage
### MATLAB
```matlab
% Generate PPM bounds table
table1
% Generate RGS plan parameters table
table23
% Compare sampling plans efficiency
table4
% Generate OC curves for n=50
fig1
% Generate OC curves for n=100
fig2
```

## File Structure
- `fig1.m`: Creates OC curves comparing single sampling and RGS plans (n=50)
- `fig2.m`: Creates OC curves comparing single sampling and RGS plans (n=100)
- `table1.m`: Generates PPM bounds for different Cpk values
- `table23.m`: Generates optimal parameters for RGS plans
- `table4.m`: Compares efficiency between single sampling and RGS plans
- `analysis1-4.m`: Additional analysis scripts for parameter studies

## Requirements
### MATLAB
- Statistics and Machine Learning Toolbox (for chi2cdf, normpdf functions)
- Optimization Toolbox (for fmincon, fsolve functions)
- MATLAB R2019b or newer recommended

## Mathematical Components
The implementation includes several key mathematical components:
- Numerical integration using F_cpk computation
- Chi-square and normal distribution calculations
- Optimization for finding optimal parameters (n, ka, kr)
- ASN calculations for efficiency comparison
- PPM bounds computation

## Implementation Notes
- Process characteristic parameter (ξ) is fixed at 1.0
- Integration tolerances: RelTol = 1e-6, AbsTol = 1e-9
- Optimization tolerances: FunctionTolerance = 1e-8, StepTolerance = 1e-8
- Maximum iterations for optimization set to 1000-2000
- Sample sizes are rounded to ensure integer values

## License
This project is available for academic and educational purposes.

## References
Wu, C.-W., Aslam, M., Chen, J.C., & Jun, C.-H. (2015). A repetitive group sampling plan by variables inspection for product acceptance determination. European Journal of Industrial Engineering, 9(3), 308-326.
