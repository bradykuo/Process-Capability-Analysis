function main_contour_plot()
    % Define parameters
    global b1 b2 xi alpha beta
    alpha = 0.05;  % Producer's risk
    beta = 0.1;    % Consumer's risk
    C_AQL = 1.33;  % Acceptable Quality Level
    C_LTPD = 1.00; % Lot Tolerance Percent Defective
    xi = 1.0;      % Using ξ = 1.0 as per paper
    b1 = 3 * C_AQL + abs(xi);
    b2 = 3 * C_LTPD + abs(xi);

    % Generate grid points
    [n_grid, c0_grid] = meshgrid(linspace(10, 150, 75), linspace(0.8, 1.5, 75));
    S1_values = zeros(size(n_grid));
    S2_values = zeros(size(n_grid));

    % Calculate S1 and S2 values with progress indicator
    fprintf('Calculating grid values...\n');
    for i = 1:size(n_grid, 1)
        if mod(i, 10) == 0
            fprintf('Processing row %d of %d\n', i, size(n_grid, 1));
        end
        for j = 1:size(n_grid, 2)
            S1_values(i,j) = calculate_S1(n_grid(i,j), c0_grid(i,j));
            S2_values(i,j) = calculate_S2(n_grid(i,j), c0_grid(i,j));
        end
    end

    % Find intersection point
    try
        [n_int, c0_int] = find_intersection();
        fprintf('Successfully found intersection point\n');
    catch ME
        fprintf('Failed to find intersection: %s\n', ME.message);
        n_int = NaN;
        c0_int = NaN;
    end

    % Create plots
    plot_results(n_grid, c0_grid, S1_values, S2_values, n_int, c0_int);
end

function result = calculate_S1(n, c0)
    global b1 xi alpha
    
    try
        upper_limit = b1 * sqrt(n);
        result = integral(@(t) S1_integrand(t, n, c0), 0, upper_limit);
        result = result - (1 - alpha);
    catch ME
        fprintf('Error in calculate_S1 for n=%f, c0=%f: %s\n', n, c0, ME.message);
        result = NaN;
    end
end

function y = S1_integrand(t, n, c0)
    global b1 xi
    G_term = chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
    phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
    y = G_term .* phi_term;
end

function result = calculate_S2(n, c0)
    global b2 xi beta
    
    try
        upper_limit = b2 * sqrt(n);
        result = integral(@(t) S2_integrand(t, n, c0), 0, upper_limit);
        result = result - beta;
    catch ME
        fprintf('Error in calculate_S2 for n=%f, c0=%f: %s\n', n, c0, ME.message);
        result = NaN;
    end
end

function y = S2_integrand(t, n, c0)
    global b2 xi
    G_term = chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
    phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
    y = G_term .* phi_term;
end

function c0_val = find_c0_for_S1(n)
    try
        fun = @(c0) calculate_S1(n, c0);
        % Test endpoints first
        f_left = fun(0.8);
        f_right = fun(1.5);
        
        if isnan(f_left) || isnan(f_right)
            fprintf('Warning: NaN values at endpoints for S1 at n=%f\n', n);
            c0_val = NaN;
            return;
        end
        
        if sign(f_left) == sign(f_right)
            fprintf('Warning: No zero crossing found for S1 at n=%f\n', n);
            c0_val = NaN;
            return;
        end
        
        c0_val = fzero(fun, [0.8, 1.5]);
    catch ME
        fprintf('Error in find_c0_for_S1 for n=%f: %s\n', n, ME.message);
        c0_val = NaN;
    end
end

function c0_val = find_c0_for_S2(n)
    try
        fun = @(c0) calculate_S2(n, c0);
        % Test endpoints first
        f_left = fun(0.8);
        f_right = fun(1.5);
        
        if isnan(f_left) || isnan(f_right)
            fprintf('Warning: NaN values at endpoints for S2 at n=%f\n', n);
            c0_val = NaN;
            return;
        end
        
        if sign(f_left) == sign(f_right)
            fprintf('Warning: No zero crossing found for S2 at n=%f\n', n);
            c0_val = NaN;
            return;
        end
        
        c0_val = fzero(fun, [0.8, 1.5]);
    catch ME
        fprintf('Error in find_c0_for_S2 for n=%f: %s\n', n, ME.message);
        c0_val = NaN;
    end
end

function [n_int, c0_int] = find_intersection()
    % First, test some points to find valid range
    n_test = linspace(10, 150, 15);
    valid_points = zeros(length(n_test), 2);
    
    for i = 1:length(n_test)
        c0_s1 = find_c0_for_S1(n_test(i));
        c0_s2 = find_c0_for_S2(n_test(i));
        valid_points(i,:) = [c0_s1, c0_s2];
        fprintf('Testing n=%f: S1 c0=%f, S2 c0=%f\n', n_test(i), c0_s1, c0_s2);
    end
    
    % Find where difference changes sign
    diffs = valid_points(:,1) - valid_points(:,2);
    valid_idx = ~isnan(diffs);
    if ~any(valid_idx)
        error('No valid intersection found in tested range');
    end
    
    % Find approximate intersection
    diffs = diffs(valid_idx);
    n_valid = n_test(valid_idx);
    sign_changes = diffs(1:end-1) .* diffs(2:end) <= 0;
    if ~any(sign_changes)
        error('No sign change found in valid range');
    end
    
    % Get index of first sign change
    idx = find(sign_changes, 1, 'first');
    n_lower = n_valid(idx);
    n_upper = n_valid(idx + 1);
    
    % Now use fzero in this narrower range
    fun = @(n) find_c0_for_S1(ceil(n)) - find_c0_for_S2(ceil(n));
    n_int = fzero(fun, [n_lower, n_upper]);
    n_int = ceil(n_int);
    c0_int = find_c0_for_S1(n_int);
end

function plot_results(n_grid, c0_grid, S1_values, S2_values, n_int, c0_int)
    figure('Position', [100 100 1200 800]);

    % Plot S1 contours
    subplot(2,2,1);
    [C,h] = contour(n_grid, c0_grid, S1_values, 10, 'b-');
    clabel(C,h,'Color','blue');
    hold on;
    contour(n_grid, c0_grid, S1_values, [0 0], 'b-', 'LineWidth', 2);
    title('Contour Plot of S1(n, c0)');
    xlabel('n');
    ylabel('c0');
    grid on;
    axis([10 150 0.8 1.5]);

    % Plot S2 contours
    subplot(2,2,2);
    [C,h] = contour(n_grid, c0_grid, S2_values, 10, 'r-');
    clabel(C,h,'Color','red');
    hold on;
    contour(n_grid, c0_grid, S2_values, [0 0], 'r-', 'LineWidth', 2);
    title('Contour Plot of S2(n, c0)');
    xlabel('n');
    ylabel('c0');
    grid on;
    axis([10 150 0.8 1.5]);

    % Combined plot
    subplot(2,1,2);
    [C1,h1] = contour(n_grid, c0_grid, S1_values, 10, 'b-');
    clabel(C1,h1,'Color','blue');
    hold on;
    [C2,h2] = contour(n_grid, c0_grid, S2_values, 10, 'r--');
    clabel(C2,h2,'Color','red');
    contour(n_grid, c0_grid, S1_values, [0 0], 'b-', 'LineWidth', 2);
    contour(n_grid, c0_grid, S2_values, [0 0], 'r-', 'LineWidth', 2);

    if ~isnan(n_int) && ~isnan(c0_int)
        plot(n_int, c0_int, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        text(n_int+5, c0_int+0.05, sprintf('(n, c0) = (%.0f, %.4f)', n_int, c0_int));
    end

    title('Combined Contour Plot of S1(n, c0) and S2(n, c0)');
    xlabel('n');
    ylabel('c0');
    legend('S1 contours', 'S2 contours', 'Location', 'best');
    grid on;
    axis([10 150 0.8 1.5]);

    % Adjust figure appearance
    set(gcf, 'Color', 'white');
    set(findall(gcf,'-property','FontSize'),'FontSize', 12);
end

main_contour_plot();
