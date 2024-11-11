% Define parameters
alpha = 0.05;  % Producer's risk
beta = 0.1;    % Consumer's risk
C_AQL = 1.33;  % Acceptable Quality Level
C_LTPD = 1.00; % Lot Tolerance Percent Defective
xi = 1.0;      % Using ξ = 1.0 as per paper
b1 = 3 * C_AQL + abs(xi);
b2 = 3 * C_LTPD + abs(xi);

% Generate grid points
n_values = linspace(10, 150, 50);
c0_values = linspace(0.8, 1.5, 50);
[N, C0] = meshgrid(n_values, c0_values);

% Initialize matrices for S1 and S2
S1 = zeros(size(N));
S2 = zeros(size(N));

% Function to calculate S1
function result = calculate_S1(n, c0, xi, b1, alpha)
    % Define integrand for S1
    function y = integrand_S1(t)
        G_term = chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
        phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = G_term .* phi_term;
    end
    
    upper_limit = b1 * sqrt(n);
    try
        result = integral(@integrand_S1, 0, upper_limit) - (1 - alpha);
    catch
        result = NaN;
    end
end

% Function to calculate S2
function result = calculate_S2(n, c0, xi, b2, beta)
    % Define integrand for S2
    function y = integrand_S2(t)
        G_term = chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
        phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = G_term .* phi_term;
    end
    
    upper_limit = b2 * sqrt(n);
    try
        result = integral(@integrand_S2, 0, upper_limit) - beta;
    catch
        result = NaN;
    end
end

% Calculate S1 and S2 for each point
for i = 1:size(N, 1)
    for j = 1:size(N, 2)
        S1(i,j) = calculate_S1(N(i,j), C0(i,j), xi, b1, alpha);
        S2(i,j) = calculate_S2(N(i,j), C0(i,j), xi, b2, beta);
    end
end

% Create figures
% Figure 1: S1 surface
figure;
surf(C0, N, S1);
colormap(winter);
title('Surface Plot of S1(n, c0)');
xlabel('c0');
ylabel('n');
zlabel('S1');
grid on;
view(45, 30);

% Figure 2: S2 surface
figure;
surf(C0, N, S2);
colormap(autumn);
title('Surface Plot of S2(n, c0)');
xlabel('c0');
ylabel('n');
zlabel('S2');
grid on;
view(45, 30);

% Figure 3: Combined surface plot
figure;
surf(C0, N, S1, 'FaceAlpha', 0.7);
hold on;
surf(C0, N, S2, 'FaceAlpha', 0.7);
% Add zero plane
surf(C0, N, zeros(size(S1)), 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3);
title('Surface Plot of S1 and S2');
xlabel('c0');
ylabel('n');
zlabel('Value');
grid on;
view(45, 30);
colormap([winter(128); autumn(128)]);
legend('S1', 'S2', 'Zero Plane');
hold off;