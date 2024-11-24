% Define process parameters
USL = 61;  % Upper Specification Limit
LSL = 40;  % Lower Specification Limit
T = 49;    % Target value

% Define different combinations of mean and standard deviation
mu = [50 52 50 52 50 52];
sigma = [2.0 2.0 3.0 3.0 3.7 3.7];

% Initialize arrays for indices
Cp = zeros(size(mu));
Cpk = zeros(size(mu));
Cpm = zeros(size(mu));

% Calculate indices for each combination
for i = 1:length(mu)
    % Calculate Cp
    Cp(i) = (USL - LSL)/(6*sigma(i));
    
    % Calculate Cpk
    Cpu = (USL - mu(i))/(3*sigma(i));
    Cpl = (mu(i) - LSL)/(3*sigma(i));
    Cpk(i) = min(Cpu, Cpl);
    
    % Calculate Cpm
    sigma_squared = sigma(i)^2 + (mu(i) - T)^2;
    Cpm(i) = (USL - LSL)/(6*sqrt(sigma_squared));
end

% Create the table
data = [mu' sigma' Cp' Cpk' Cpm'];
columnNames = {'μ', 'σ', 'Cp', 'Cpk', 'Cpm'};
T = array2table(data, 'VariableNames', columnNames);

% Display the table
disp(T)