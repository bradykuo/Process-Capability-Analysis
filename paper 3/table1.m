% Function to calculate PPM bounds for given Cpk
function [lower_ppm, upper_ppm] = calculate_ppm_bounds(cpk)
    % For one-sided specification:
    upper_ppm = (1 - normcdf(3 * cpk)) * 1e6;
    
    % For two-sided specification:
    % Lower bound is when process is centered (symmetrical case)
    lower_ppm = 2 * (1 - normcdf(3 * cpk)) * 1e6;
end

% Define Cpk values
cpk_values = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.24, 1.25, ...
              1.30, 1.33, 1.40, 1.45, 1.50, 1.60, 1.67, 1.70, 1.80, 1.90, 2.00];

% Calculate PPM bounds for each Cpk value
results = zeros(length(cpk_values), 2);
for i = 1:length(cpk_values)
    [results(i,1), results(i,2)] = calculate_ppm_bounds(cpk_values(i));
end

% Round the results to 3 decimal places
results = round(results, 3);

% Create and display the table
T = table(cpk_values', results(:,1), results(:,2), ...
    'VariableNames', {'Cpk', 'Lower_Bound', 'Upper_Bound'});
disp(T)

% Verify some values against the paper
fprintf('\nVerification of selected values:\n');
[lower_1_33, upper_1_33] = calculate_ppm_bounds(1.33);
fprintf('For Cpk = 1.33: Expected ≈ 66/33 PPM, Calculated: %.3f %.3f PPM\n', ...
    lower_1_33, upper_1_33);

[lower_1_67, upper_1_67] = calculate_ppm_bounds(1.67);
fprintf('For Cpk = 1.67: Expected ≈ 0.544/0.272 PPM, Calculated: %.3f %.3f PPM\n', ...
    lower_1_67, upper_1_67);