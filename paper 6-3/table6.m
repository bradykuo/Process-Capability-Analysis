function table6()
    % Parameters
    alphas = [0.01, 0.05, 0.10];
    betas = [0.01, 0.05, 0.10];
    quality_levels = [
        1.33, 1.00;
        1.50, 1.00;
        1.67, 1.33;
        2.00, 1.67
    ];
    m = 3; % Fixed number of stages
    
    % Store results
    results = [];
    
    % Loop through all combinations
    for i = 1:size(quality_levels, 1)
        cAQL = quality_levels(i, 1);
        cLQL = quality_levels(i, 2);
        for alpha = alphas
            for beta = betas
                % Optimize parameters
                [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m);
                % Store results
                results = [results; alpha, beta, cAQL, cLQL, n, ka, kr];
            end
        end
    end
    
    % Display results in table format
    T = array2table(results, 'VariableNames', {'alpha', 'beta', 'cAQL', 'cLQL', 'n', 'ka', 'kr'});
    disp(T);
end

function [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m)
    % Initial guess based on paper values
    x0 = [30, 1.3, 1.0]; % [n, ka, kr]
    
    % Bounds
    lb = [20, 1.0, 0.8];
    ub = [500, 3.0, 3.0];
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxIterations', 2000);
    
    % Objective function
    function obj = objective(x)
        n_val = x(1);
        ka_val = x(2);
        kr_val = x(3);
        ASN_AQL = calc_ASN(cAQL, n_val, ka_val, kr_val, m, 0);  % xi=0 for AQL
        ASN_LQL = calc_ASN(cLQL, n_val, ka_val, kr_val, m, 1);  % xi=1 for LQL
        obj = 0.5 * (ASN_AQL + ASN_LQL);
    end
    
    % Constraints
    function [c, ceq] = constraints(x)
        n_val = x(1);
        ka_val = x(2);
        kr_val = x(3);
        piA_AQL = calc_piA(cAQL, n_val, ka_val, kr_val, m, 0);  % xi=0 for AQL
        piA_LQL = calc_piA(cLQL, n_val, ka_val, kr_val, m, 1);  % xi=1 for LQL
        c = [(1-alpha) - piA_AQL;   % piA(cAQL) >= 1-alpha
             piA_LQL - beta;        % piA(cLQL) <= beta
             kr_val - ka_val];      % ka > kr
        ceq = [];
    end
    
    % Run optimization
    [x_opt, ~] = fmincon(@objective, x0, [], [], [], [], lb, ub, @constraints, options);
    
    % Return optimized parameters
    n = ceil(x_opt(1));
    ka = x_opt(2);
    kr = x_opt(3);
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