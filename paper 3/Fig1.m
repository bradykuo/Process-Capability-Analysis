% Function to calculate PPM bounds for given Cpk
function [lower_ppm, upper_ppm] = calculate_ppm_bounds(cpk)
    % For one-sided specification:
    upper_ppm = (1 - normcdf(3 * cpk)) * 1e6;
    
    % For two-sided specification:
    % Lower bound is when process is centered (symmetrical case)
    lower_ppm = 2 * (1 - normcdf(3 * cpk)) * 1e6;
end

% Create sequence of Cpk values from 0.6 to 2.0
cpk_seq = 0.6:0.01:2.0;

% Calculate PPM bounds for each Cpk value
results = zeros(length(cpk_seq), 2);
for i = 1:length(cpk_seq)
    [results(i,1), results(i,2)] = calculate_ppm_bounds(cpk_seq(i));
end

% Create the plot
figure;
plot(cpk_seq, results(:,1), 'LineWidth', 1);
hold on;
plot(cpk_seq, results(:,2), 'LineWidth', 1);

% Set axis labels and limits
xlabel('C_{pk}', 'FontSize', 12);
ylabel('PPM', 'FontSize', 12);
ylim([0 80000]);
xlim([0.6 2.0]);

% Set custom y-axis ticks
yticks(0:20000:80000);

% Add subtitle
text(0.5, -15000, 'Fig. 1. The bounds on nonconforming units in PPM versus C_{pk}.', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% Adjust plot margins
set(gca, 'Position', [0.15 0.2 0.8 0.7]);

hold off;