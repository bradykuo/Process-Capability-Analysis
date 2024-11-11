% Function to calculate the integral in equations (9) and (10)
function result = calculate_integral(n, c0, C_value, xi)
    b = 3 * C_value + abs(xi);
    
    % Define integrand function
    function y = integrand(t)
        % G is the chi-square CDF with n-1 degrees of freedom
        G_term = chi2cdf((n - 1) * (b * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
        % phi is the standard normal PDF
        phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = G_term .* phi_term;
    end
    
    % Numerical integration
    result = integral(@integrand, 0, b*sqrt(n));
end

% Function to find c0 for given n and parameters
function [c0_opt, fval] = find_c0(n, C_AQL, C_LTPD, alpha, beta, xi)
    % Objective function
    function err = objective(c0)
        eq1 = calculate_integral(n, c0, C_AQL, xi) - (1 - alpha);
        eq2 = calculate_integral(n, c0, C_LTPD, xi) - beta;
        err = eq1^2 + eq2^2;
    end
    
    % Find optimal c0
    [c0_opt, fval] = fminbnd(@objective, 0.8, 2.0);
end

% Main script
% Generate data for plots
xi_values = 0:0.1:2;
CAQL_values = [1.33, 1.50, 1.67, 2.00];
CLTPD = 1.00;
alpha = 0.05;
beta = 0.05;

% Initialize arrays
n_results = zeros(length(xi_values), length(CAQL_values));
c0_results = zeros(length(xi_values), length(CAQL_values));

% Calculate values for each CAQL and xi combination
for i = 1:length(CAQL_values)
    CAQL = CAQL_values(i);
    for j = 1:length(xi_values)
        xi = xi_values(j);
        
        % Find minimum n that satisfies both equations
        n_min = 10;  % Start with small n
        while true
            [c0, fval] = find_c0(n_min, CAQL, CLTPD, alpha, beta, xi);
            
            % Check if equations are satisfied
            eq1 = calculate_integral(n_min, c0, CAQL, xi) >= (1 - alpha);
            eq2 = calculate_integral(n_min, c0, CLTPD, xi) <= beta;
            
            if eq1 && eq2
                break
            end
            n_min = n_min + 1;
            if n_min > 100
                break  % Safety break
            end
        end
        
        % Store results
        n_results(j,i) = n_min;
        c0_results(j,i) = c0;
    end
end

% Create plots
figure('Position', [100 100 1200 500]);

% Plot (a) - Required Sample Size
subplot(1,2,1);
plot(xi_values, n_results, 'LineWidth', 2);
ylim([0 100]);
grid on;
title('Required Sample Size vs \xi');
xlabel('\xi');
ylabel('Required sample size n');
legend(cellstr(num2str(CAQL_values', 'C_{AQL}=%.2f')), 'Location', 'best');

% Plot (b) - Critical Acceptance Value
subplot(1,2,2);
plot(xi_values, c0_results, 'LineWidth', 2);
ylim([1.0 1.6]);
grid on;
title('Critical Acceptance Value vs \xi');
xlabel('\xi');
ylabel('Critical value C_o');
legend(cellstr(num2str(CAQL_values', 'C_{AQL}=%.2f')), 'Location', 'best');

% Adjust figure appearance
set(gcf, 'Color', 'white');
set(findall(gcf,'-property','FontSize'),'FontSize', 12);