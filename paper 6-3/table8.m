function table8()
    % Parameters for first case (α,β)=(0.01,0.05)
    alpha1 = 0.01;
    beta1 = 0.05;
    
    % Parameters for second case (α,β)=(0.05,0.10)
    alpha2 = 0.05;
    beta2 = 0.10;
    
    % Common parameters
    cAQL = 1.33;
    cLQL = 1.00;
    
    % Generate Cpk values
    Cpk = 0.80:0.05:1.50;
    
    % Initialize results matrix
    results = zeros(length(Cpk), 9);  % 9 columns: Cpk, SSP1, m2_1, m3_1, m4_1, SSP2, m2_2, m3_2, m4_2
    results(:,1) = Cpk';
    
    % Calculate for first case (α,β)=(0.01,0.05)
    [n_ssp1, ka_ssp1, kr_ssp1] = optimize_parameters(alpha1, beta1, cAQL, cLQL, 1);
    [n2_1, ka2_1, kr2_1] = optimize_parameters(alpha1, beta1, cAQL, cLQL, 2);
    [n3_1, ka3_1, kr3_1] = optimize_parameters(alpha1, beta1, cAQL, cLQL, 3);
    [n4_1, ka4_1, kr4_1] = optimize_parameters(alpha1, beta1, cAQL, cLQL, 4);
    
    % Calculate for second case (α,β)=(0.05,0.10)
    [n_ssp2, ka_ssp2, kr_ssp2] = optimize_parameters(alpha2, beta2, cAQL, cLQL, 1);
    [n2_2, ka2_2, kr2_2] = optimize_parameters(alpha2, beta2, cAQL, cLQL, 2);
    [n3_2, ka3_2, kr3_2] = optimize_parameters(alpha2, beta2, cAQL, cLQL, 3);
    [n4_2, ka4_2, kr4_2] = optimize_parameters(alpha2, beta2, cAQL, cLQL, 4);
    
    % Calculate ASN values
    for i = 1:length(Cpk)
        % For (α,β)=(0.01,0.05)
        results(i,2) = n_ssp1;  % SSP ASN is just n
        results(i,3) = calc_ASN(Cpk(i), n2_1, ka2_1, kr2_1, 2, 1);
        results(i,4) = calc_ASN(Cpk(i), n3_1, ka3_1, kr3_1, 3, 1);
        results(i,5) = calc_ASN(Cpk(i), n4_1, ka4_1, kr4_1, 4, 1);
        
        % For (α,β)=(0.05,0.10)
        results(i,6) = n_ssp2;  % SSP ASN is just n
        results(i,7) = calc_ASN(Cpk(i), n2_2, ka2_2, kr2_2, 2, 1);
        results(i,8) = calc_ASN(Cpk(i), n3_2, ka3_2, kr3_2, 3, 1);
        results(i,9) = calc_ASN(Cpk(i), n4_2, ka4_2, kr4_2, 4, 1);
    end
    
    % Display results
    fprintf('\nα      β     Cpk    SSP     m=2     m=3     m=4\n');
    fprintf('------------------------------------------------\n');
    
    % First group (α,β)=(0.01,0.05)
    for i = 1:length(Cpk)
        fprintf('0.01   0.05   %.2f   %3d   %6.2f  %6.2f  %6.2f\n', ...
            results(i,1), results(i,2), results(i,3), results(i,4), results(i,5));
    end
    
    % Second group (α,β)=(0.05,0.10)
    for i = 1:length(Cpk)
        fprintf('0.05   0.10   %.2f   %3d   %6.2f  %6.2f  %6.2f\n', ...
            results(i,1), results(i,6), results(i,7), results(i,8), results(i,9));
    end
end

function [n, ka, kr] = optimize_parameters(alpha, beta, cAQL, cLQL, m)
    % Multiple initial guesses
    initial_points = [
        [100, 1.3, 1.0];
        [150, 1.5, 1.3];
        [80, 1.2, 1.1];
        [120, 1.4, 1.2]
    ];
    
    % Bounds
    lb = [20, 1.0, 0.8];
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
    
    % Try each initial point
    best_obj = Inf;
    best_x = [];
    
    for i = 1:size(initial_points, 1)
        [x_tmp, fval] = fmincon(@objective, initial_points(i,:), [], [], [], [], lb, ub, @constraints, options);
        
        if fval < best_obj
            best_obj = fval;
            best_x = x_tmp;
        end
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