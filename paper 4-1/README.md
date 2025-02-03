# Process Capability Testing Implementation
Implementation of statistical testing procedures for one-sided specification limits based on process capability indices CPU and CPL, as described in [Lin & Pearn (2002)](https://doi.org/10.1016/S0026-2714(02)00103-8).

## Features
- Calculates p-values and critical values for process capability testing
- Generates probability density function (PDF) plots for different CPU values
- Computes non-conforming units in parts per million (NCPPM)
- Supports multiple sample sizes and significance levels
- Handles one-sided specification limits (both upper and lower)
- Provides comprehensive reference tables for practical implementation

## Installation
```bash
git clone https://github.com/yourusername/process-capability-testing
cd process-capability-testing
```

## Usage
### MATLAB
```matlab
% Generate Table 1 (NCPPM values)
table1
% Generate Table 2 (Critical values for C=1.25)
table2
% Generate Table 3 (Critical values for C=1.45)
table3
% Generate Table 4 (Critical values for C=1.60)
table4
% Generate PDF plots
fig1to4
```

## File Structure
- `table1.m`: Generates NCPPM values for different CPU values
- `table2.m`: Generates critical values for C=1.25
- `table3.m`: Generates critical values for C=1.45
- `table4.m`: Generates critical values for C=1.60
- `fig1to4.m`: Creates PDF plots for different CPU values
- `appendixA.m`: Calculates p-values for hypothesis testing
- `appendixB.m`: Calculates critical values for hypothesis testing

## Requirements
### MATLAB
- Statistics and Machine Learning Toolbox (for normcdf, nctcdf, nctinv functions)
- MATLAB R2019b or newer recommended

## Mathematical Components
The implementation includes several key mathematical components:
- Non-central t-distribution calculations
- Numerical integration for PDF computation
- Gamma function calculations for correction factors
- Normal distribution probability computations

## Implementation Notes
- Uses log-gamma function for numerical stability
- Implements error handling for computational edge cases
- Provides formatted output matching paper tables
- Sample sizes range from 10 to 505 in steps of 5

## License
This project is available for academic and educational purposes.

## References
Lin, P. C., & Pearn, W. L. (2002). Testing process capability for one-sided specification limit with application to the voltage level translator. Microelectronics Reliability, 42(12), 1975-1983.
