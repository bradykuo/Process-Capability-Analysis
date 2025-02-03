# Variables Stage-Independent Multiple Sampling Plan Implementation

Implementation of a Stage-Independent Multiple Sampling Plan (SIMSP) based on the process capability index Cpk for lot determination. This project provides MATLAB implementations for optimizing sampling parameters, generating Operating Characteristic (OC) curves, and calculating Average Sample Numbers (ASN).

## Features
- Optimizes sampling plan parameters (n, ka, kr) for multiple stages (m=1 to m=10)
- Generates OC curves comparing SSP and SIMSP with different m values
- Calculates and plots ASN curves for various quality levels
- Computes actual producer's and consumer's risks (α*, β*)
- Supports multiple quality level combinations (cAQL/cLQL)
- Provides comprehensive parameter tables for practical implementation
- Includes robustness checks and numerical stability improvements

## Installation
```bash
git clone https://github.com/yourusername/simsp-implementation
cd simsp-implementation
```

## Usage
### MATLAB
```matlab
% Generate ASN curves for first quality level combination
fig1ab
fig1cd

% Generate ASN curves for second quality level combination
fig2ab
fig2cd

% Generate OC curves
fig3a
fig3b

% Generate parameter tables
table1  % Parameters for (α,β)=(0.01,0.05)
table2  % Parameters for (α,β)=(0.05,0.05)
table3  % Parameters for (α,β)=(0.05,0.10)
table4  % Parameters for (α,β)=(0.10,0.10)

% Generate risk analysis tables
table10  % Actual risks for (α,β)=(0.01,0.05)
table11  % Actual risks for (α,β)=(0.05,0.10)
```

## File Structure
- `fig1ab.m`, `fig1cd.m`: ASN curves for (cAQL,cLQL)=(1.33,1.00)
- `fig2ab.m`, `fig2cd.m`: ASN curves for (cAQL,cLQL)=(2.00,1.67)
- `fig3a.m`, `fig3b.m`: OC curves for different risk combinations
- `table1.m` to `table4.m`: Parameter tables for various α,β combinations
- `table5.m` to `table7.m`: Plan parameters for m=2,3,4
- `table8.m`, `table9.m`: ASN comparison tables
- `table10.m`, `table11.m`: Risk analysis tables
- `table13.m`: Process yield analysis

## Requirements
### MATLAB
- Optimization Toolbox (for fmincon function)
- Statistics and Machine Learning Toolbox (for normcdf, normpdf functions)
- MATLAB R2019b or newer recommended

## Mathematical Components
- CDF calculation for Cpk estimator (F_Cpk function)
- Probability of resampling (PS) calculation
- Operating characteristic function
- Average sample number computation
- Chi-square and normal distribution integration
- Sequential quadratic programming optimization

## Implementation Notes
- The F_Cpk function uses numerical integration with RelTol=1e-6 and AbsTol=1e-9
- Multiple initial points are used for robust optimization
- Constraints include producer's risk, consumer's risk, and acceptance criteria
- Edge cases are handled in probability calculations for numerical stability
- Sample sizes are rounded up to ensure integer values
- Optimization uses both 'sqp' and 'interior-point' algorithms depending on the context

## License
This project is available for academic and educational purposes.

## References
Wu, C.-W., Darmawan, A., & Liu, S.-W. (2022). Stage-independent multiple sampling plan by variables inspection for lot determination based on the process capability index Cpk. International Journal of Production Research. 
