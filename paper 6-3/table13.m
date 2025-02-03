function table13()
    % Fixed parameters
    alpha = 0.05;
    beta = 0.10;
    cAQL = 1.33;
    cLQL = 1.00;
    
    % Cpk values as given in the table
    Cpk = [0.80; 0.90; 1.00; 1.14; 1.18; 1.20; 1.24; 1.28; 1.30; 1.34];
    
    % Get optimized parameters for m=3 (SIMSP)
    [n_ssp, ka_ssp, kr_ssp] = optimize_parameters(alpha, beta, cAQL, cLQL, 1);
    [n_simsp, ka_simsp, kr_simsp] = optimize_parameters(alpha, beta, cAQL, cLQL, 3);
    
    % Calculate process yield (lower bound)
    process_yield = 2 * normcdf(3*Cpk) - 1;
    process_yield = process_yield * 100;  % Convert to percentage
    
    % Calculate ASN for SIMSP for each Cpk value
    ASN_values = zeros(size(Cpk));
    for i = 1:length(Cpk)
        ASN_values(i) = calc_ASN(Cpk(i), n_simsp, ka_simsp, kr_simsp, 3, 1);
    end
    
    % Display header
    fprintf('\nTable 13. The required ASN for Cpk-based SSP and SIMSP under\n');
    fprintf('the given conditions, (cAQL, cLQL) = (1.33, 1.00) and (α, β) = (0.05, 0.10).\n\n');
    
    % Column headers
    fprintf('%-6s %-15s %-14s %-20s\n', 'Cpk', 'Process yield', 'Cpk-based SSP', 'The proposed Cpk-based SIMSP');
    fprintf('%-6s %-15s %-14s %-10s %-10s\n', '', '(lower bound)', 'n', 'n', 'ASN');
    fprintf('%s\n', repmat('-', 1, 65));
    
    % Print results
    for i = 1:length(Cpk)
        fprintf('%.2f   %8.4f%%      %3d           %3d          %6.2f\n', ...
            Cpk(i), process_yield(i), n_ssp, n_simsp, ASN_values(i));
    end
end

function [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m)
    % Initial points based on expected solution ranges
    initial_points = [
        [175, 1.85, 1.75];
        [185, 1.90, 1.80];
        [195, 1.95, 1.85];
        [165, 1.80, 1.70];
        [155, 1.75, 1.65]
    ];
    
    % Set bounds based on problem context
    lb = [30, 1.0, 0.9];
    ub = [300, 2.5, 2.4];
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point', ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'OptimalityTolerance', 1e-8, ...
        'MaxFunctionEvaluations', 2000, ...
        'MaxIterations', 1000, ...
        'ConstraintTolerance', 1e-6);
    
    % Try each initial point
    best_obj = Inf;
    best_x = [];
    best_feas = Inf;
    
    for i = 1:size(initial_points, 1)
        try
            [x_tmp, fval] = fmincon(@objective, initial_points(i,:), [], [], [], [], lb, ub, @constraints, options);
            
            % Check constraints
            [c, ~] = constraints(x_tmp);
            feas = max(max(c, 0));
            
            % Update best solution if more feasible or equally feasible but better objective
            if feas < best_feas || (abs(feas - best_feas) < 1e-6 && fval < best_obj)
                best_feas = feas;
                best_obj = fval;
                best_x = x_tmp;
            end
        catch
            continue;
        end
    end
    
    if isempty(best_x)
        best_x = initial_points(1,:);  % Use first initial point as fallback
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
        
        c = [(1-alpha) - piA_AQL;    % Producer's risk constraint
             piA_LQL - beta;         % Consumer's risk constraint
             kr_val - ka_val + 0.1;  % Ensure kr < ka with minimum gap
             ka_val - kr_val - 0.3]; % Maximum gap between ka and kr
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
    if abs(PS - 1) < 1e-10
        piA = Pa;
    else
        piA = Pa * (1 - PS^m) / (1 - PS);
    end
end

function ASN = calc_ASN(c, n, ka, kr, m, xi)
    PS = calc_PS(c, n, ka, kr, xi);
    if abs(PS - 1) < 1e-10
        ASN = n;
    else
        ASN = n * (1 - PS^m) / (1 - PS);
    end
end