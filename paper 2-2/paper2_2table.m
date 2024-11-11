function calculatePanels()
    % Start timing
    tic;
    
    % Define parameters
    n_values = 10:5:200;
    c_pk_hat_values_A = 0.7:0.1:1.8;
    c_pk_hat_values_B = 1.9:0.1:3.0;
    
    % Calculate panels
    fprintf('Calculating Panel A...\n');
    panel_A = zeros(length(n_values), length(c_pk_hat_values_A));
    for i = 1:length(c_pk_hat_values_A)
        for j = 1:length(n_values)
            panel_A(j,i) = round(calculate_lcb(n_values(j), c_pk_hat_values_A(i)), 3);
        end
    end
    
    fprintf('Calculating Panel B...\n');
    panel_B = zeros(length(n_values), length(c_pk_hat_values_B));
    for i = 1:length(c_pk_hat_values_B)
        for j = 1:length(n_values)
            panel_B(j,i) = round(calculate_lcb(n_values(j), c_pk_hat_values_B(i)), 3);
        end
    end
    
    % Create tables
    result_table_A = array2table([n_values' panel_A], 'VariableNames', ...
        ['n', cellstr(num2str(c_pk_hat_values_A', '%.1f'))']);
    result_table_B = array2table([n_values' panel_B], 'VariableNames', ...
        ['n', cellstr(num2str(c_pk_hat_values_B', '%.1f'))']);
    
    % Display results
    disp('Panel A:');
    disp(result_table_A);
    disp('Panel B:');
    disp(result_table_B);
    
    % Display execution time
    execution_time = toc;
    fprintf('\nTotal execution time: %.2f seconds\n', execution_time);
end

function lcb = calculate_lcb(n, Cpk_hat, xi_hat, gamma)
    if nargin < 3
        xi_hat = 1.0;
    end
    if nargin < 4
        gamma = 0.95;
    end
    
    % Define options for fzero
    options = optimset('Display', 'off');
    
    % Find root using fzero
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

calculatePanels()