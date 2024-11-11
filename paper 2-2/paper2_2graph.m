function plotCpkAnalysis()
    % Start timing
    tic;
    
    % Parameters
    n_values = [30, 50, 70, 100, 150, 200];
    xi_values = 0:0.1:3;
    Cpk_hat_values = [0.7, 0.9, 1.2, 2.0, 2.5, 3.0];
    
    % Create figure
    figure('Position', [100, 100, 1200, 800]);
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Color map for different n values
    colors = lines(length(n_values));
    
    % Create plots
    for i = 1:length(Cpk_hat_values)
        results = calculate_C_vs_xi(Cpk_hat_values(i), n_values, xi_values);
        
        nexttile
        hold on
        for j = 1:length(n_values)
            idx = results.n == n_values(j);
            plot(results.xi(idx), results.C(idx), 'Color', colors(j,:), 'LineWidth', 1.5)
        end
        hold off
        
        title(['Ĉpk = ', num2str(Cpk_hat_values(i))])
        xlabel('|ξ|')
        ylabel('C')
        grid on
        box on
    end
    
    % Add common title
    title(t, 'C vs |ξ| for different Ĉpk values', 'FontSize', 16, 'FontWeight', 'bold')
    
    % Add common legend
    leg = legend(string(n_values), 'Location', 'southoutside', 'Orientation', 'horizontal');
    leg.Layout.Tile = 'south';
    title(leg, 'n')
    
    % Display execution time
    execution_time = toc;
    fprintf('Total execution time: %.2f seconds\n', execution_time);
end

function results = calculate_C_vs_xi(Cpk_hat, n_values, xi_values)
    % Initialize results structure
    [N, Xi] = meshgrid(n_values, xi_values);
    results.n = N(:);
    results.xi = Xi(:);
    results.C = zeros(size(results.n));
    
    % Calculate C for each combination
    for i = 1:length(results.n)
        results.C(i) = calculate_lcb(results.n(i), Cpk_hat, results.xi(i));
    end
end

function lcb = calculate_lcb(n, Cpk_hat, xi_hat, gamma)
    if nargin < 4
        gamma = 0.95;
    end
    
    % Use fzero to find the root
    options = optimset('Display', 'off');
    lcb = fzero(@(C) integrate_func(C, n, Cpk_hat, xi_hat, gamma), [0, Cpk_hat], options);
end

function result = integrate_func(C, n, Cpk_hat, xi_hat, gamma)
    b = 3*C + abs(xi_hat);
    
    % Define integrand function
    integrand = @(t) G_func(t, n, b, Cpk_hat) .* phi_func(t, xi_hat, n);
    
    % Perform numerical integration
    result = integral(integrand, 0, b*sqrt(n)) - (1 - gamma);
end

function G = G_func(t, n, b, Cpk_hat)
    G = chi2cdf((n-1)*(b*sqrt(n)-t).^2 ./ (9*n*Cpk_hat^2), n-1);
end

function phi = phi_func(t, xi_hat, n)
    phi = normpdf(t + xi_hat*sqrt(n)) + normpdf(t - xi_hat*sqrt(n));
end

plotCpkAnalysis()