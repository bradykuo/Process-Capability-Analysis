# Bootstrap Confidence Limits for Process Capability Indices

## Overview
This project implements bootstrap methods to calculate confidence limits for process capability indices (Cp, Cpk, and Cpm) under different probability distributions. The implementation is available in both MATLAB and R, providing flexibility for different user preferences.

## Features
- Calculates three types of bootstrap confidence limits:
  - Standard Bootstrap (SB)
  - Percentile Bootstrap (PB)
  - Bias-Corrected Percentile Bootstrap (BCPB)
- Supports multiple probability distributions:
  - Normal distribution (table1, table2)
  - Log-normal distribution (table3)
  - Chi-squared distribution (table4)
- Provides coverage probability analysis
- Implements parallel simulations for efficiency

## Requirements

### MATLAB Implementation
- MATLAB (R2019b or later recommended)
- Statistics and Machine Learning Toolbox

### R Implementation
- R (version 3.6.0 or later recommended)
- Base R packages (no additional packages required)

## File Structure
```
.
├── matlab/
│   ├── table1.m
│   ├── table2.m
│   ├── table3.m
│   └── table4.m
└── R/
    ├── table1.R
    ├── table2.R
    ├── table3.R
    └── table4.R
```

## Usage

### MATLAB
```matlab
% Run normal distribution simulation
table1
% Run bootstrap analysis with normal distribution
table2
% Run lognormal distribution simulation
table3
% Run chi-squared distribution simulation
table4
```

### R
```r
# Run normal distribution simulation
source("table1.R")
# Run bootstrap analysis with normal distribution
source("table2.R")
# Run lognormal distribution simulation
source("table3.R")
# Run chi-squared distribution simulation
source("table4.R")
```

## Implementation Details

### Process Capability Indices
1. **Cp (Process Capability)**
   ```
   Cp = (USL - LSL)/(6*σ)
   ```

2. **Cpk (Process Capability with mean shift)**
   ```
   Cpk = min((USL - μ)/(3*σ), (μ - LSL)/(3*σ))
   ```

3. **Cpm (Taguchi Capability Index)**
   ```
   Cpm = (USL - LSL)/(6*sqrt(σ² + (μ - T)²))
   ```

### Bootstrap Methods

#### Standard Bootstrap (SB)
- Uses normal approximation
- Calculates bootstrap standard deviation
- Confidence limit: estimate ± Z_α * bootstrap_std

#### Percentile Bootstrap (PB)
- Uses ordered bootstrap replicates
- Takes percentiles directly
- More robust to non-normality

#### Bias-Corrected Percentile Bootstrap (BCPB)
- Accounts for potential bias
- Uses normal score transformation
- Adjusts percentiles based on observed bias

## Parameters
- USL (Upper Specification Limit) = 61
- LSL (Lower Specification Limit) = 40
- T (Target Value) = 49
- α (Significance Level) = 0.05
- B (Bootstrap Resamples) = 1,000
- N (Simulation Repetitions) = 1,000

## Output Format
The programs output coverage probabilities for each bootstrap method:
```
mu: XX.X, sigma: XX.X, n: XX
True Values - Cp: XX.XX, Cpk: XX.XX, Cpm: XX.XX
Coverage Probabilities:
        Cp      Cpk     Cpm
SB:     X.XXX   X.XXX   X.XXX
PB:     X.XXX   X.XXX   X.XXX
BCPB:   X.XXX   X.XXX   X.XXX
```

## Performance Considerations
- Large number of bootstrap resamples (B=1000) ensures stability
- Multiple sample sizes (n=20, 40, 70) for comprehensive analysis
- Different distribution parameters for thorough testing

## License
This project is available for academic and educational purposes.

## References
1. Efron, B. (1979). Bootstrap methods: another look for the jack-knife. The Annals of Statistics
2. Franklin, L.A. and Wasserman, G.S. (1992). Bootstrap lower confidence limits for process capability indices
3. Efron, B. and Tibshirani, R. (1986). Bootstrap methods for standard errors, confidence intervals

## Acknowledgments
Special thanks to the original paper authors and the statistical community for developing these methods.
