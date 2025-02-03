% Clear workspace and command window
clear;
clc;

% Input parameters
C = 1.25;  % Quality requirement
n = 120;   % Sample size
w = 1.433; % Estimated CPU

% Calculate degrees of freedom and non-centrality parameter
F = n - 1;
ND = 3 * sqrt(n) * C;

% Calculate DN (same as SAS program)
DN = sqrt((n-2)/2) * (1 - 1/(4*(n-2)) + ...
    1/(32*(n-2)^2) + 5/(128*(n-2)^3));

% Calculate BF
BF = sqrt(2/(n-1)) * DN;

% Calculate X
X = 3 * sqrt(n) * w/BF;

% Calculate P-value (using non-central t-distribution)
PV = 1 - nctcdf(X, F, ND);

% Create table for output
T = table(C, n, w, PV, 'VariableNames', {'C', 'n', 'w', 'PV'});

% Display results
disp('The output is:')
disp('MATLAB System')
disp(' ')
disp(T)

% Format output to match SAS exactly
fprintf('\nFormatted output:\n');
fprintf('Obs    C        n        w        PV\n');
fprintf('1      %.2f    %d     %.3f    %.3f\n', C, n, w, PV);