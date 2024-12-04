% Initialize parameters
alpha_values = [0.01, 0.025, 0.05, 0.075, 0.10];  % Producer's risk α
beta_values = [0.01, 0.025, 0.05, 0.075, 0.10];   % Consumer's risk β

% Set CAQL and CLTPD parameter 
parameter_sets = [
    1.33, 1.00;   
    1.50, 1.00;   
    1.50, 1.33;   
    1.67, 1.33;   
    1.67, 1.50;  
    2.00, 1.67    
];

xi = 0;  % Distribution characteristic parameter (according to paper)

% Initialize results matrix
results = [];

% Loop through all combinations
for k = 1:size(parameter_sets, 1)
    CAQL = parameter_sets(k, 1);
    CLTPD = parameter_sets(k, 2);
    
    for i = 1:length(alpha_values)
        alpha = alpha_values(i);
        for j = 1:length(beta_values)
            beta = beta_values(j);
            
            % Calculate b1 and b2 according to paper's formula
            b1 = 3 * CAQL * (1 + xi^2)^0.5;
            b2 = 3 * CLTPD * (1 + xi^2)^0.5;
            
            % Define S1 and S2 functions
            S1 = @(n, c0) integral(@(t) chi2cdf((b1.^2.*n)./(9.*c0.^2) - t.^2, n-1) .* ...
                (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
                0, b1.*sqrt(n)./(3.*c0)) - (1-alpha);
            
            S2 = @(n, c0) integral(@(t) chi2cdf((b2.^2.*n)./(9.*c0.^2) - t.^2, n-1) .* ...
                (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
                0, b2.*sqrt(n)./(3.*c0)) - beta;
            
            % Define system of equations
            system_of_equations = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
            
            % Set initial guesses scaled by CAQL
            initial_n = 50 * CAQL/1.33;      
            initial_c0 = 1.1 * CAQL/1.33;    
            initial_guess = [initial_n, initial_c0];
            
            % Solver options
            options = optimoptions('fsolve', ...
                'Display', 'off', ...
                'FunctionTolerance', 1e-10, ...
                'StepTolerance', 1e-10, ...
                'MaxIterations', 1000);
            
            % Solve equations
            try
                solution = fsolve(system_of_equations, initial_guess, options);
                n = ceil(solution(1));
                c0 = solution(2);
                
                % Validate solution
                if n > 0 && n < 1000 && c0 > 1.0 && c0 < 2.0
                    results = [results; alpha, beta, CAQL, CLTPD, n, c0];
                else
                    results = [results; alpha, beta, CAQL, CLTPD, NaN, NaN];
                end
            catch
                results = [results; alpha, beta, CAQL, CLTPD, NaN, NaN];
            end
            
            % Display progress
            fprintf('α=%.3f, β=%.3f, CAQL=%.2f, CLTPD=%.2f: n=%d, c0=%.4f\n', ...
                alpha, beta, CAQL, CLTPD, n, c0);
        end
    end
end

% Convert results to table format
result_table = array2table(results, ...
    'VariableNames', {'Alpha', 'Beta', 'CAQL', 'CLTPD', 'n', 'c0'});

% Display results
disp(result_table);