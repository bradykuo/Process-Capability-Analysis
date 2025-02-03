% Define CPU values
CPU = [0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 2.1 2.3];

% Calculate PPM values
% For CPU, use right tail probability (1-normcdf)
% multiply by 1e6 to convert to PPM
PPM = (1 - normcdf(3*CPU))*1e6;

% Display results in table format
fprintf('Various CPU values and the corresponding nonconformities\n');
fprintf('----------------------------------------------------\n');
fprintf('CPU\t\tppm\n');
fprintf('----------------------------------------------------\n');
for i = 1:length(CPU)
    fprintf('%.1f\t\t%.7f\n', CPU(i), PPM(i));
end