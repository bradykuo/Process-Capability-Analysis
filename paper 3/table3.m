function results = generate_table3()
    % Define parameters
    alpha_values = [0.01, 0.025, 0.05, 0.075, 0.10];
    beta_values = [0.01, 0.025, 0.05, 0.075, 0.10];
    cases = {
        struct('c_aql', 1.33, 'c_ltpd', 1.00),
        struct('c_aql', 1.50, 'c_ltpd', 1.33),
        struct('c_aql', 1.67, 'c_ltpd', 1.33),
        struct('c_aql', 2.00, 'c_ltpd', 1.67)
    };
    
    % Create results matrix
    num_rows = length(alpha_values) * length(beta_values);
    results = zeros(num_rows, 10);  % alpha, beta + 4 pairs of (n, c0)
    
    row = 1;
    for alpha = alpha_values
        for beta = beta_values
            results(row, 1:2) = [alpha, beta];
            
            for i = 1:length(cases)
                case_params = cases{i};
                [n, c0] = find_n_c0(alpha, beta, case_params.c_aql, case_params.c_ltpd);
                results(row, (2*i+1):(2*i+2)) = [n, c0];
            end
            
            row = row + 1;
            fprintf('Calculated: alpha = %.3f, beta = %.3f\n', alpha, beta);
        end
    end
    
    % Print formatted table
    fprintf('\nTable 3: Required sample sizes (n) and critical acceptance values (c0)\n');
    fprintf('α      β       CAQL=1.33      CAQL=1.50      CAQL=1.67      CAQL=2.00\n');
    fprintf('               CLTPD=1.00      CLTPD=1.33      CLTPD=1.33      CLTPD=1.67\n');
    fprintf('               n      c0       n      c0       n      c0       n      c0\n');
    fprintf('%s\n', repmat('-', 1, 75));
    
    for i = 1:size(results, 1)
        fprintf('%.3f  %.3f  ', results(i,1), results(i,2));
        for j = 1:4
            fprintf('%3d  %.4f  ', results(i,2*j+1), results(i,2*j+2));
        end
        fprintf('\n');
    end
end

function [best_n, best_c0] = find_n_c0(alpha, beta, c_aql, c_ltpd)
    % Set search range
    n_min = 30;
    n_max = 1000;
    
    best_n = NaN;
    best_c0 = NaN;
    min_error = Inf;
    
    for n = n_min:n_max
        c0 = find_c0_for_n(n, alpha, beta, c_aql, c_ltpd);
        error = evaluate_n(n, alpha, beta, c_aql, c_ltpd);
        
        if error < min_error
            min_error = error;
            best_n = n;
            best_c0 = c0;
        end
    end
    
    best_c0 = round(best_c0 * 10000) / 10000;  % Round to 4 decimal places
end

function c0 = find_c0_for_n(n, alpha, beta, c_aql, c_ltpd)
    % Define objective function
    objective = @(c0) abs(calculate_S1(n, c0, c_aql, alpha)) + ...
                     abs(calculate_S2(n, c0, c_ltpd, beta));
    
    % Use fminbnd instead of optimize
    [c0, ~] = fminbnd(objective, c_ltpd, c_aql + 0.5);
end

function error = evaluate_n(n, alpha, beta, c_aql, c_ltpd)
    c0 = find_c0_for_n(n, alpha, beta, c_aql, c_ltpd);
    s1 = calculate_S1(n, c0, c_aql, alpha);
    s2 = calculate_S2(n, c0, c_ltpd, beta);
    error = abs(s1) + abs(s2);
end

function result = calculate_S1(n, c0, c_aql, alpha, xi)
    if nargin < 5
        xi = 1;
    end
    
    b1 = 3 * c_aql + abs(xi);
    
    % Define integrand function
    function val = integrand(t)
        G_val = chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 ./ (9 * n * c0^2), n - 1);
        phi_val = normpdf(t + sqrt(n)) + normpdf(t - sqrt(n));
        val = G_val .* phi_val;
    end
    
    % Numerical integration using integral
    result = integral(@integrand, 0, b1 * sqrt(n)) - (1 - alpha);
end

function result = calculate_S2(n, c0, c_ltpd, beta, xi)
    if nargin < 5
        xi = 1;
    end
    
    b2 = 3 * c_ltpd + abs(xi);
    
    % Define integrand function
    function val = integrand(t)
        G_val = chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 ./ (9 * n * c0^2), n - 1);
        phi_val = normpdf(t + sqrt(n)) + normpdf(t - sqrt(n));
        val = G_val .* phi_val;
    end
    
    % Numerical integration using integral
    result = integral(@integrand, 0, b2 * sqrt(n)) - beta;
end