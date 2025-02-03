% Clear workspace and command window
clear;
clc;

% Input parameters
C = 1.25;      % Capability requirement
Alpha = 0.05;  % Risk level
n = 120;       % Sample size

% Calculate degrees of freedom and non-centrality parameter
F = n - 1;
ND = 3 * sqrt(n) * C;

% Calculate DN (same as SAS program)
DN = sqrt((n-2)/2) * (1 - 1/(4*(n-2)) + ...
    1/(32*(n-2)^2) + 5/(128*(n-2)^3));

% Calculate BF
BF = sqrt(2/(n-1)) * DN;

% Calculate Critical Value c0
% Note: In MATLAB, we use nctinv instead of TINV for non-central t distribution
C0 = BF/(3 * sqrt(n)) * nctinv(1-Alpha, F, ND);

% Create table for output
T = table(C, Alpha, n, C0, 'VariableNames', {'C', 'Alpha', 'n', 'C0'});

% Display results
disp('The output is:')
disp('MATLAB System')
disp(' ')
disp(T)

% Format output to match SAS exactly
fprintf('\nFormatted output:\n');
fprintf('Obs    C        Alpha     n        C0\n');
fprintf('1      %.2f    %.3f     %d      %.3f\n', C, Alpha, n, C0);
