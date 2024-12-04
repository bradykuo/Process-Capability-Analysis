# Variables Sampling Plan Implementation

This repository implements a variables sampling plan based on the process capability index Cpm for quality control applications. The implementation provides tools for determining optimal sample sizes and critical acceptance values while ensuring protection for both producers and consumers.

## Overview

The sampling plan helps determine:
- Required sample size (n) for inspection
- Critical acceptance value (C₀) for lot sentencing
- Operating characteristics that satisfy both producer's and consumer's risks

## Key Features

- Calculation of optimal sampling parameters
- Visualization tools including contour and surface plots
- Comprehensive lookup tables for various quality levels
- Support for different risk levels (α and β)
- Consideration of process capability indices

## Files Description

- `contour_plot.m`: Generates contour plots showing acceptance regions and optimal parameter combinations
- `figure5.m`: Creates 3D surface plots demonstrating relationships between risks and sampling parameters
- `surface_plot.m`: Visualizes acceptance surfaces and their intersections
- `table1.m`: Generates lookup tables for various parameter combinations

## Mathematical Background

The implementation is based on the following key equations:

1. Process Capability Index (Cpm):
```
Cpm = (USL - LSL) / (6 * sqrt(σ² + (μ - T)²))
```
where:
- USL, LSL: Upper and lower specification limits
- μ: Process mean
- T: Target value
- σ: Process standard deviation

2. Acceptance Criteria:
```
b₁ = 3 * CAQL * (1 + ξ²)^0.5
b₂ = 3 * CLTPD * (1 + ξ²)^0.5
```

where ξ is the distribution characteristic parameter.

## Usage

1. Set your quality parameters:
```matlab
alpha = 0.05;     % Producer's risk
beta = 0.05;      % Consumer's risk
CAQL = 1.33;      % Acceptable Quality Level
CLTPD = 1.00;     % Lot Tolerance Percent Defective
```

2. Run the desired analysis script:
```matlab
% For contour plots
contour_plot

% For surface visualization
surface_plot

% For lookup tables
table1
```

## Installation Requirements

- MATLAB R2018b or later
- Optimization Toolbox
- Statistics and Machine Learning Toolbox

## Input Parameters

- `alpha`: Producer's risk (Type I error)
- `beta`: Consumer's risk (Type II error)
- `CAQL`: Acceptable Quality Level
- `CLTPD`: Lot Tolerance Percent Defective
- `xi`: Distribution characteristic parameter

## Output

1. Numerical Results:
- Optimal sample size (n)
- Critical acceptance value (C₀)

2. Visualizations:
- Contour plots of acceptance regions
- 3D surface plots of risk relationships
- Combined visualization of acceptance criteria

## Example Results

For typical values (α = 0.05, β = 0.05, CAQL = 1.33, CLTPD = 1.00):
- Required sample size (n) ≈ 68
- Critical acceptance value (C₀) ≈ 1.1668

## Applications

This implementation is suitable for:
- Quality control in manufacturing
- Process capability analysis
- Acceptance sampling procedures
- Quality improvement initiatives

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is available for academic and educational purposes.

## References

1. Pearn, W.L., & Wu, C.W. (2006). Variables sampling plans with PPM fraction of defectives and process loss consideration. Journal of the Operational Research Society, 57, 450-459.
