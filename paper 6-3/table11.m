function table11()
    % Fixed parameters
    alpha = 0.05;
    beta = 0.10;
    
    % Quality level combinations
    quality_levels = [
        1.33, 1.00;
        1.50, 1.00;
        1.67, 1.33;
        2.00, 1.67
    ];
    
    % For each m value (2,3,4)
    fprintf('\nm    (cAQL,cLQL)      n      ka      kr       α*       β*\n');
    fprintf('-----------------------------------------------------------\n');
    
    for m = 2:4
        for i = 1:size(quality_levels,1)
            cAQL = quality_levels(i,1);
            cLQL = quality_levels(i,2);
            
            % Get n, ka, kr using optimize_parameters
            [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m);
            
            % Calculate α* (using ξ=0 for AQL)
            alpha_star = 1 - calc_piA(cAQL, n, ka, kr, m, 0);
            
            % Calculate β* (using ξ=1 for LQL)
            beta_star = calc_piA(cLQL, n, ka, kr, m, 1);
            
            % Print results
            fprintf('%d    (%.2f,%.2f)     %3d   %.4f  %.4f  %.5f  %.5f\n', ...
                m, cAQL, cLQL, n, ka, kr, alpha_star, beta_star);
        end
    end
end

function [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m)
    % More initial guesses with wider spread
    initial_points = [
        [100, 1.3, 1.0];    % Original
        [150, 1.5, 1.3];
        [80, 1.2, 1.1];
        [120, 1.4, 1.2];
        [200, 1.6, 1.4];    % Higher values
        [180, 1.7, 1.5];
        [160, 1.8, 1.6];
        [140, 1.9, 1.7];
        [90, 1.3, 1.2];     % Lower values
        [70, 1.4, 1.3];
        [60, 1.5, 1.4];
        [50, 1.6, 1.5]
    ];
    
    % Adjusted bounds
    lb = [30, 1.0, 0.9];    % Increased lower bounds
    ub = [300, 2.5, 2.4];   % Decreased upper bounds
    
    % Improved optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxFunctionEvaluations', 20000, ...
        'MaxIterations', 20000, ...
        'ScaleProblem', 'obj-and-constr');
    
    % Try each initial point
    best_obj = Inf;
    best_x = [];
    best_feas = Inf;  % Track constraint violation
    
    for i = 1:size(initial_points, 1)
        try
            [x_tmp, fval, exitflag, output, lambda] = fmincon(@objective, initial_points(i,:), [], [], [], [], lb, ub, @constraints, options);
            
            % Calculate constraint violation
            [c, ~] = constraints(x_tmp);
            feas = max(max(c, 0));
            
            % Update best solution if more feasible or equally feasible but better objective
            if feas < best_feas || (feas == best_feas && fval < best_obj)
                best_feas = feas;
                best_obj = fval;
                best_x = x_tmp;
            end
        catch
            continue;
        end
    end
    
    % Check if we found a solution
    if isempty(best_x)
        error('No feasible solution found');
    end
    
    % Return optimized parameters
    n = ceil(best_x(1));
    ka = best_x(2);
    kr = best_x(3);
    
    function obj = objective(x)
        n_val = x(1);
        ka_val = x(2);
        kr_val = x(3);
        ASN_AQL = calc_ASN(cAQL, n_val, ka_val, kr_val, m, 0);
        ASN_LQL = calc_ASN(cLQL, n_val, ka_val, kr_val, m, 1);
        obj = 0.5 * (ASN_AQL + ASN_LQL);
    end
    
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