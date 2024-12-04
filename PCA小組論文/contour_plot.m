% Define ranges
n_values = 10:1:150;          
c0_values = 0.8:0.01:1.5;     

% Initialize matrices for S1 and S2 values
S1_values = zeros(length(c0_values), length(n_values));
S2_values = zeros(length(c0_values), length(n_values));

% Define parameters
alpha = 0.05;     % Producer's risk
beta = 0.05;      % Consumer's risk
CAQL = 1.33;      % Acceptable Quality Level
CLTPD = 1.0;      % Lot Tolerance Percent Defective
xi = 0;           % Distribution characteristic parameter

% Calculate b1 and b2 according to equations
b1 = 3*CAQL*(1 + xi^2)^0.5;
b2 = 3*CLTPD*(1 + xi^2)^0.5;

% Define S1 and S2 functions (outside loop for efficiency)
S1 = @(n, c0) integral(@(t) chi2cdf((b1.^2.*n)./(9.*c0.^2) - t.^2, n-1) ...
    .* (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
    0, b1.*sqrt(n)./(3.*c0), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1-alpha);

S2 = @(n, c0) integral(@(t) chi2cdf((b2.^2.*n)./(9.*c0.^2) - t.^2, n-1) ...
    .* (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
    0, b2.*sqrt(n)./(3.*c0), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;

% Calculate S1 and S2 values for each combination
for i = 1:length(c0_values)
    c0 = c0_values(i);
    for j = 1:length(n_values)
        n = n_values(j);
        S1_values(i, j) = S1(n, c0);
        S2_values(i, j) = S2(n, c0);
    end
    % Display progress
    if mod(i, 10) == 0
        fprintf('Processing... %.1f%%\n', i/length(c0_values)*100);
    end
end

[N, C0] = meshgrid(n_values, c0_values);

% Plot Figure 1: S1 Contour Plot
figure;
contour_levels_S1 = -0.9:0.1:0;
contour(N, C0, S1_values, contour_levels_S1, 'ShowText', 'on');
xlabel('Sample Size (n)');
ylabel('Critical Value (C_0)');
title('S_1(n, C_0) Contour Plot');
axis([10 150 0.8 1.5]);
box on;
grid off;

% Plot Figure 2: S2 Contour Plot
figure;
contour_levels_S2 = 0:0.1:0.8;
contour(N, C0, S2_values, contour_levels_S2, 'ShowText', 'on');
xlabel('Sample Size (n)');
ylabel('Critical Value (C_0)');
title('S_2(n, C_0) Contour Plot');
axis([10 150 0.8 1.5]);
box on;
grid off;

% Define objective function for solving
obj_fun = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];

% Multiple initial guesses for better convergence
initial_guesses = [
    68, 1.1668;   
    50, 1.2;      
    60, 1.15      
];

% Find best solution from multiple initial guesses
best_solution = [];
best_fval_norm = Inf;

for i = 1:size(initial_guesses, 1)
    try
        [x_sol, fval] = fsolve(obj_fun, initial_guesses(i,:), ...
            optimoptions('fsolve', 'Display', 'off', ...
            'FunctionTolerance', 1e-10, ...
            'StepTolerance', 1e-10));
        
        current_fval_norm = norm(fval);
        if current_fval_norm < best_fval_norm
            best_fval_norm = current_fval_norm;
            best_solution = x_sol;
        end
    catch
        continue;
    end
end

% Process solution results
if ~isempty(best_solution)
    n_intersect = best_solution(1);
    c0_intersect = best_solution(2);
    fprintf('Found solution: n = %.4f, C0 = %.4f\n', n_intersect, c0_intersect);
else
    warning('No valid solution found');
    n_intersect = NaN;
    c0_intersect = NaN;
end

% Plot Figure 3: Combined Contour Plot with Intersection Point
figure;
hold on;

% Plot all contour lines (dashed)
[C1,h1] = contour(N, C0, S1_values, contour_levels_S1, '--k');
[C2,h2] = contour(N, C0, S2_values, contour_levels_S2, '--k');

% Plot special lines for S1=0 and S2=0
[C3,h3] = contour(N, C0, S1_values, [0 0], 'b-', 'LineWidth', 2);
[C4,h4] = contour(N, C0, S2_values, [0 0], 'r-', 'LineWidth', 2);

% Plot intersection point if found
if ~isnan(n_intersect)
    plot(n_intersect, c0_intersect, 'k.', 'MarkerSize', 15);
    text(n_intersect+5, c0_intersect, sprintf('(n, C_0) = (%.0f, %.4f)', ...
        n_intersect, c0_intersect), 'FontSize', 10);
end

xlabel('Sample Size (n)');
ylabel('Critical Value (C_0)');
title('Contour Plot with Intersection Point');
axis([10 150 0.8 1.5]);
box on;
grid off;
hold off;