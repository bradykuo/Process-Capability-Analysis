% Define ranges 
n_values = 150:-1:10;        
c0_values = 1.5:-0.01:0.8;    

% Initialize matrices
S1_values = zeros(length(c0_values), length(n_values));
S2_values = zeros(length(c0_values), length(n_values));

% Define parameters
alpha = 0.05;     % Producer's risk
beta = 0.05;      % Consumer's risk
CAQL = 1.33;      % Acceptable Quality Level
CLTPD = 1.0;      % Lot Tolerance Percent Defective
xi = 0;           % Distribution characteristic parameter

% Calculate b1 and b2
b1 = 3 * CAQL + xi;
b2 = 3 * CLTPD + xi;

% Calculate S1 and S2 for each (n, c0) combination
for i = 1:length(c0_values)
    c0 = c0_values(i);
    for j = 1:length(n_values)
        n = n_values(j);
        
        % S1 integral
        S1 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 / (9 * n * c0^2), n-1) ...
            .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b1 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1 - alpha);
        
        % S2 integral
        S2 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 / (9 * n * c0^2), n-1) ...
            .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b2 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;
        
        S1_values(i, j) = S1(n, c0);
        S2_values(i, j) = S2(n, c0);
    end
end

[N, C0] = meshgrid(n_values, c0_values);

% Fig 1: S1 surface plot
figure;
surf(N, C0, S1_values, 'EdgeColor', 'none');
xlabel('Sample Size (n)');
ylabel('Critical Value (c0)');
zlabel('S1(n, c0)');
title('Surface Plot of S1(n, c0)');
colorbar;
view(45, 30);    % Adjust viewing angle for better surface visualization

% Fig 2: S2 surface plot
figure;
surf(N, C0, S2_values, 'EdgeColor', 'none');
xlabel('Sample Size (n)');
ylabel('Critical Value (c0)');
zlabel('S2(n, c0)');
title('Surface Plot of S2(n, c0)');
colorbar;
view(45, 30);    % Adjust viewing angle for better surface visualization

% Fig 3: Combined surface plot
figure;
hold on;
surf(N, C0, S1_values, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'b');
surf(N, C0, S2_values, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'r');
xlabel('Sample Size (n)');
ylabel('Critical Value (c0)');
zlabel('S1 and S2 values');
title('3D Surface Plot of S1(n, c0) and S2(n, c0)');
legend({'S1(n, c0)', 'S2(n, c0)'}, 'Location', 'best');
colorbar;
view(45, 30);    % Adjust viewing angle for better surface visualization
hold off;