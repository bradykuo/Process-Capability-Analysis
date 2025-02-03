function table4()
    % Parameters to compute table for
    alphas = [0.01, 0.05, 0.1];
    betas = [0.01, 0.05, 0.1];
    quality_levels = {[1.33, 1.00], [1.5, 1.33], [1.67, 1.33], [2.00, 1.67]};
    
    % Create arrays to store results
    results = [];
    
    for ql = 1:length(quality_levels)
        CAQL = quality_levels{ql}(1);
        CLTPD = quality_levels{ql}(2);
        for alpha = alphas
            for beta = betas
                % Calculate single sampling plan
                xi = 1.0;
                % 計算規格界限係數
                b1 = 3 * CAQL + xi;
                b2 = 3 * CLTPD + xi;
                % 定義積分函數
                S1 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
                    .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b1 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1 - alpha);
                S2 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
                    .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b2 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;
                % 建立聯立方程組
                system_of_equations = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
                % 設定初始值和求解選項
                initial_guess = [50 + 20 * (CAQL - 1), 1.0 + 0.5 * (CAQL - 1)];
                options = optimoptions('fsolve', ...
                    'Display', 'off', ...
                    'FunctionTolerance', 1e-8, ...
                    'StepTolerance', 1e-8, ...
                    'MaxIterations', 1000);
                % 求解方程組
                solution = fsolve(system_of_equations, initial_guess, options);
                n_single = ceil(solution(1)); % 樣本數取整數
            
                % Calculate RGS plan
                x = find_optimal_parameters(alpha, beta, CAQL, CLTPD);
                n = x(1);
                ka = x(2);
                kr = x(3);
                
                % Calculate ASN using CLTPD
                bL = 3*CLTPD + abs(xi);
                asn = n / (Pa(n, ka, bL) + Pr(n, kr, bL));
                
                % Store both results
                results = [results; CAQL, CLTPD, alpha, beta, n_single, round(asn)];
            end
        end
    end
    
    % Display results
    headers = {'CAQL', 'CLTPD', 'α', 'β', 'Variables single sampling plan', 'Variables RGS plan'};
    disp(array2table(results, 'VariableNames', headers));
end

function x = find_optimal_parameters(alpha, beta, CAQL, CLTPD)
    % Fixed starting point
    x0 = [30, 1.3, 1.0];
    
    % Define bounds
    lb = [20, 1.0, 0.8];
    ub = [500, 3.0, 3.0];
    
    % Define optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxIterations', 2000);
    
    % Calculate b values
    xi = 1.0;
    bA = 3*CAQL + abs(xi);
    bL = 3*CLTPD + abs(xi);
    
    % Define minimum ASN objective function
    min_ASN = @(x) x(1) / (Pa(x(1), x(2), bL) + Pr(x(1), x(3), bL));
    
    % Define nonlinear constraints
    constraints = @(x) deal([], [
        (1 - alpha) - (Pa(x(1), x(2), bA) / (Pa(x(1), x(2), bA) + Pr(x(1), x(3), bA)));
        (Pa(x(1), x(2), bL) / (Pa(x(1), x(2), bL) + Pr(x(1), x(3), bL))) - beta
    ]);
    
    % Run optimization
    x = fmincon(min_ASN, x0, [], [], [], [], lb, ub, constraints, options);
end

function pa = Pa(n, ka, b)
    xi = 1.0;
    pa = compute_F_cpk(ka, n, b, xi);
end

function pr = Pr(n, kr, b)
    xi = 1.0;
    pr = 1 - compute_F_cpk(kr, n, b, xi);
end

function F = compute_F_cpk(ka, n, b, xi)
    % Define integrand function
    function y = integrand(t)
        chi_term = (n-1) * (b*sqrt(n) - t).^2 / (9*n*ka^2);
        phi_terms = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = chi2cdf(chi_term, n-1) .* phi_terms;
    end
    
    % Compute integral
    F = integral(@integrand, 0, b*sqrt(n));
end