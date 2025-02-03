% Define Cpk index values as shown in table
cpk_values = [0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.24 1.25 1.30 1.33 1.40 1.45 1.50 1.60 1.67 1.70 1.80 1.90 2.00];

% Function to calculate PPM bounds
function [upper_ppm, lower_ppm] = calculatePPMBounds(Cpk)
    % Calculate yields
    upper_yield = 100 * (2 * normcdf(3*Cpk) - 1);
    lower_yield = 100 * normcdf(3*Cpk);
    
    % Convert to PPM
    upper_ppm = (100 - lower_yield) * 10000;
    lower_ppm = (100 - upper_yield) * 10000;
end

% Calculate and display results in table format
fprintf('Index            Lower bound           Upper bound      \n');
fprintf('----------------------------------------------------------\n');

for i = 1:length(cpk_values)
    [lower_ppm, upper_ppm] = calculatePPMBounds(cpk_values(i));
    
    % Format numbers according to table style
    if lower_ppm >= 14
        lower_str = sprintf('%8.0f', lower_ppm);  
    else
        lower_str = sprintf('%8.3f', lower_ppm);  
    end
    
    if upper_ppm >= 13
        upper_str = sprintf('%8.0f', upper_ppm);  
    else
        upper_str = sprintf('%8.3f', upper_ppm);  
    end
    
    fprintf('%5.2f      %15s      %15s\n', cpk_values(i), lower_str, upper_str);
end