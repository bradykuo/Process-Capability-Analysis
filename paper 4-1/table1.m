% Define CPU values
CPU = [0.80; 1.00; 1.25; 1.45; 1.60; 1.80; 2.00];

% Calculate NCPPM values
% For CPU, NCPPM = (1 - normcdf(3*CPU)) * 1e6
NCPPM = (1 - normcdf(3*CPU)) * 1e6;

% Create table
T = table(CPU, NCPPM, 'VariableNames', {'CPU', 'NCPPM'});

% Display table with formatted output
disp('CPU and the corresponding non-conforming units (in ppm)')
disp('---------------------------------------------------')
fprintf('Values of CPU       NCPPM\n')
disp('---------------------------------------------------')
for i = 1:height(T)
    if T.NCPPM(i) >= 1
        fprintf('%6.2f            %.0f\n', T.CPU(i), T.NCPPM(i))
    else
        fprintf('%6.2f            %.5f\n', T.CPU(i), T.NCPPM(i))
    end
end
disp('---------------------------------------------------')
