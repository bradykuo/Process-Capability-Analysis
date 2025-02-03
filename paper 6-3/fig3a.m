function fig3a()
    % Parameters
    alpha = 0.01;
    beta = 0.05;
    cAQL = 1.33;
    cLQL = 1.00;
    
    % Get optimized parameters for different m values
    [n_ssp, ka_ssp, kr_ssp] = optimize_parameters(alpha, beta, cAQL, cLQL, 1);
    [n2, ka2, kr2] = optimize_parameters(alpha, beta, cAQL, cLQL, 2);
    [n3, ka3, kr3] = optimize_parameters(alpha, beta, cAQL, cLQL, 3);
    [n4, ka4, kr4] = optimize_parameters(alpha, beta, cAQL, cLQL, 4);
    
    % Generate Cpk values
    Cpk = 0.9:0.005:1.4;
    
    % Calculate OC (Probability of acceptance) for each plan
    Pa_SSP = zeros(size(Cpk));
    Pa_m2 = zeros(size(Cpk));
    Pa_m3 = zeros(size(Cpk));
    Pa_m4 = zeros(size(Cpk));
    
    % Calculate Pa values
    for i = 1:length(Cpk)
        % For SSP (m=1)
        PS_ssp = calc_PS(Cpk(i), n_ssp, ka_ssp, kr_ssp, 1);
        Pa = 1 - F_Cpk(ka_ssp, n_ssp, Cpk(i), 1);
        Pa_SSP(i) = Pa * (1 - PS_ssp) / (1 - PS_ssp);  % For m=1, same as Pa
        
        % For m=2
        PS_m2 = calc_PS(Cpk(i), n2, ka2, kr2, 1);
        Pa = 1 - F_Cpk(ka2, n2, Cpk(i), 1);
        Pa_m2(i) = Pa * (1 - PS_m2^2) / (1 - PS_m2);
        
        % For m=3
        PS_m3 = calc_PS(Cpk(i), n3, ka3, kr3, 1);
        Pa = 1 - F_Cpk(ka3, n3, Cpk(i), 1);
        Pa_m3(i) = Pa * (1 - PS_m3^3) / (1 - PS_m3);
        
        % For m=4
        PS_m4 = calc_PS(Cpk(i), n4, ka4, kr4, 1);
        Pa = 1 - F_Cpk(ka4, n4, Cpk(i), 1);
        Pa_m4(i) = Pa * (1 - PS_m4^4) / (1 - PS_m4);
    end
    
    % Create plot
    figure('Position', [100 100 800 600]);
    
    % Plot curves
    plot(Cpk, Pa_m2, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(Cpk, Pa_m3, 'b--', 'LineWidth', 1.5);
    plot(Cpk, Pa_m4, 'Color', '#006400', 'LineStyle', ':', 'LineWidth', 1.5);
    plot(Cpk, Pa_SSP, 'k-', 'LineWidth', 1.5);
    
    % Customize plot
    grid off;
    xlabel('C_{pk}', 'FontSize', 12);
    ylabel('Probability of acceptance', 'FontSize', 12);
    ylim([0 1]);
    xlim([0.9 1.4]);
    
    % Add legend
    legend('The proposed plan (m=2)', 'The proposed plan (m=3)', ...
           'The proposed plan (m=4)', 'Variables SSP', ...
           'Location', 'southeast', 'FontSize', 10);
    
    % Add title
    title('(c_{AQL}, c_{LQL}) = (1.33,1.00) and (\alpha, \beta) = (0.01, 0.05)', 'FontSize', 12);
    
    % Fine-tune appearance
    set(gca, 'Box', 'on');
    set(gca, 'LineWidth', 1);
    set(gca, 'FontSize', 10);
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