function c0 = calculate_c0_fminbnd(C, Cp, n, alpha)
    % Initial interval [a, b]
    a = 0.5;
    b = 5;
    
    % Solve using fminbnd
    [c0, fval, exitflag] = fminbnd(@(x) abs(objective_function(x, C, Cp, n, alpha)), a, b);
    
    % Check if solution was found
    if exitflag == 1
        fprintf('Calculated c0: %f\n', c0);
    else
        error('fminbnd did not converge to a solution.');
    end
end

function F = objective_function(c0, C, Cp, n, alpha)
    % Define the integrand
    integrand = @(y) G_term(y, n, Cp, c0) .* f_term(y, C, Cp, n);
    
    % Compute the integral using integral function
    try
        integral_result = integral(integrand, 0, 3*Cp*sqrt(n));
        F = integral_result - alpha;
    catch
        error('Integration failed. Check the function or parameter values.');
    end
end

function G = G_term(y, n, Cp, c0)
    % Calculate chi-square CDF term
    G = chi2cdf((n-1)*(3*Cp*sqrt(n) - y).^2 ./ (9*n*c0^2), n-1);
end

function f = f_term(y, C, Cp, n)
    % Calculate folded normal density term
    f = normpdf(y + 3*(Cp - C)*sqrt(n)) + normpdf(y - 3*(Cp - C)*sqrt(n));
end

% Example usage
C = 1.33;
n = 50;
alpha = 0.05;
Cp = C + (n < 100)*0.33 + (n >= 100)*0.12;
c0 = calculate_c0_fminbnd(C, Cp, n, alpha);