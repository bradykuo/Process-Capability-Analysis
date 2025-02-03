function table23()
    % Parameters to compute table for
    alphas = [0.01, 0.025, 0.05, 0.075, 0.1];
    betas = [0.01, 0.025, 0.05, 0.075, 0.1];
    quality_levels = {[1.33, 1.00], [1.50, 1.33], [1.67, 1.33], [2.00, 1.67]};
    
    % Create arrays to store results
    results = [];
    
    for ql = 1:length(quality_levels)
        CAQL = quality_levels{ql}(1);
        CLTPD = quality_levels{ql}(2);
        for alpha = alphas
            for beta = betas
                x = find_optimal_parameters(alpha, beta, CAQL, CLTPD);
                n = x(1);
                ka = x(2);
                kr = x(3);
                % Calculate ASN using CLTPD
                xi = 1.0;
                bL = 3*CLTPD + abs(xi);
                asn = n / (Pa(n, ka, bL) + Pr(n, kr, bL));
                results = [results; CAQL, CLTPD, alpha, beta, round(n), ka, kr, round(asn)];
            end
        end
    end
    
    % Display results
    headers = {'CAQL', 'CLTPD', 'α', 'β', 'n', 'ka', 'kr', 'ASN'};
    disp(array2table(results, 'VariableNames', headers));
end

function x = find_optimal_parameters(alpha, beta, CAQL, CLTPD)
    % Initial guess and bounds
    x0 = [30, 1.3, 1.0];
    lb = [20, 1.0, 0.8];
    ub = [500, 3.0, 3.0];
    
    % Define optimization options
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-9, ...
        'MaxIterations', 1000);

    % Calculate b values
    xi = 1.0;
    bA = 3*CAQL + abs(xi);  % For AQL
    bL = 3*CLTPD + abs(xi); % For LTPD
    
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