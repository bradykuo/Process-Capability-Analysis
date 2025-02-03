# Acceptance Sampling Plans Implementation

Implementation of an effective decision-making method for product acceptance based on process capability indices, as described in [Pearn & Wu (2007)](https://doi.org/10.1016/j.omega.2005.01.018).

## Features

- Calculates minimum sample sizes and critical values for acceptance sampling
- Generates contour and surface plots of acceptance probabilities
- Computes PPM bounds for various Cpk values
- Handles both single and double-sided specification limits

## Installation

```bash
git clone https://github.com/yourusername/process-capability-analysis
cd process-capability-analysis
```

## Usage

### MATLAB
```matlab
% Calculate minimum sample size
min_n

% Generate Table 3
table3
```

### R
```r
source("Table3.R")
table3 <- generate_table3()
```

## File Structure

- `min_n.m`: Calculates minimum sample size
- `table3.m/R`: Generates sampling plan parameters
- `contour_plot.m/R`: Creates contour plots
- `surface_plot.m/R`: Creates 3D surface plots
- `Fig1.m/R`: Generates PPM vs Cpk plots

## Requirements

### MATLAB
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

### R
- pracma
- ggplot2
- dplyr
- plotly

## License

This project is available for academic and educational purposes.

## References

Pearn, W.L., Wu, C.-W. (2007). An effective decision making method for product acceptance. Omega, 35(1), 12-21.
