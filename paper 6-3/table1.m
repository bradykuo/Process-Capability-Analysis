function table1()
    % Fixed parameters
    alpha = 0.01;
    beta = 0.05;
    quality_levels = [
        1.33, 1.00;
        1.50, 1.00;
        1.67, 1.33;
        2.00, 1.67
    ];
    
    % Initialize results array
    results = [];
    
    % Calculate for each combination
    for i = 1:size(quality_levels, 1)
        cAQL = quality_levels(i,1);
        cLQL = quality_levels(i,2);
        
        for m = 1:10
            fprintf('Calculating for cAQL=%.2f, cLQL=%.2f, m=%d\n', cAQL, cLQL, m);
            
            % Optimize parameters
            [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m);
            
            % Calculate ASNMax at midpoint Cpk = (ka+kr)/2
            Cpk_mid = (ka + kr)/2;
            ASNMax = calc_ASN(Cpk_mid, n, ka, kr, m, 1);  % Use xi=1 for conservative estimate
            
            % Store results
            results = [results; cAQL, cLQL, m, n, ASNMax];
        end
    end
    
    % Display results in table format
    fprintf('\nResults for (α,β) = (0.01,0.05)\n');
    fprintf('CAQL   CLQL    m    n     ASNMax\n');
    fprintf('----------------------------------\n');
    for i = 1:size(results,1)
        fprintf('%.2f   %.2f   %2d   %3d   %.2f\n', ...
            results(i,1), results(i,2), results(i,3), results(i,4), results(i,5));
    end
end

function [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m)
    % Multiple initial guesses
    x0_list = [
        [100, 1.3, 1.0];  % Original guess
        [80, 1.2, 1.1];   % Another guess
        [50, 1.4, 1.2];   % Third guess
        [120, 1.5, 1.3]   % Fourth guess
    ];
    
    % Bounds
    lb = [10, 1.0, 0.8];
    ub = [500, 3.0, 3.0];
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 10000);
    
    % Try multiple initial points
    best_obj = Inf;
    best_x = [];
    
    for i = 1:size(x0_list, 1)
        % Try current initial point
        [x_opt, fval, exitflag] = fmincon(@objective, x0_list(i,:), [], [], [], [], lb, ub, @constraints, options);
        
        % Check if this solution is better and feasible
        if exitflag > 0 && fval < best_obj
            [c, ~] = constraints(x_opt);
            if all(c <= 1e-6)  % Check if constraints are satisfied
                best_obj = fval;
                best_x = x_opt;
            end
        end
    end
    
    % If no solution found, throw error
    if isempty(best_x)
        error('No feasible solution found');
    end
    
    % Return optimized parameters
    n = ceil(best_x(1));
    ka = best_x(2);
    kr = best_x(3);
    
    % Nested objective function
    function obj = objective(x)
        n_val = x(1);
        ka_val = x(2);
        kr_val = x(3);
        ASN_AQL = calc_ASN(cAQL, n_val, ka_val, kr_val, m, 0);
        ASN_LQL = calc_ASN(cLQL, n_val, ka_val, kr_val, m, 1);
        obj = 0.5 * (ASN_AQL + ASN_LQL);
    end
    
    % Nested constraints function
    function [c, ceq] = constraints(x)
        n_val = x(1);
        ka_val = x(2);
        kr_val = x(3);
        piA_AQL = calc_piA(cAQL, n_val, ka_val, kr_val, m, 0);
        piA_LQL = calc_piA(cLQL, n_val, ka_val, kr_val, m, 1);
        c = [(1-alpha) - piA_AQL;
            piA_LQL - beta;
            kr_val - ka_val];
        ceq = [];
    end
end

function F = F_Cpk(y, n, c, xi)
    b = 3*c + abs(xi);
    
    function val = integrand(t)
        chi_term = (n-1)*(b*sqrt(n)-t).^2/(9*n*y^2);
        phi_terms = normpdf(t - xi*sqrt(n)) + normpdf(t + xi*sqrt(n));
        val = chi2cdf(chi_term, n-1) .* phi_terms;
    end
    
    F = 1 - integral(@integrand, 0, b*sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9);
end

function PS = calc_PS(c, n, ka, kr, xi)
    PS = F_Cpk(ka, n, c, xi) - F_Cpk(kr, n, c, xi);
end

function piA = calc_piA(c, n, ka, kr, m, xi)
    Pa = 1 - F_Cpk(ka, n, c, xi);
    PS = calc_PS(c, n, ka, kr, xi);
    piA = Pa * (1 - PS^m) / (1 - PS);
end

function ASN = calc_ASN(c, n, ka, kr, m, xi)
    PS = calc_PS(c, n, ka, kr, xi);
    ASN = n * (1 - PS^m) / (1 - PS);
end