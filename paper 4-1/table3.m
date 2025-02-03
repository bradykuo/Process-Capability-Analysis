% Clear the workspace and command window
clear all;
clc;

% Function to calculate bn-1 using the log gamma method
function b = calculate_bn1_log(n)
    % Use log gamma to avoid overflow
    log_b = 0.5*log(2/(n-1)) + gammaln((n-1)/2) - gammaln((n-2)/2);
    b = exp(log_b);
end

% Main program starts
function calculate_critical_values()
    % Set parameters
    C = 1.45;  % Change process capability requirement 
    alpha_values = [0.01, 0.025, 0.05];  % Significance levels
    n_values = 10:5:505;  % Range of sample sizes

    % Initialize results matrix
    results = zeros(length(n_values), length(alpha_values));

    % Calculate critical values
    for i = 1:length(n_values)
        n = n_values(i);
        try
            % Calculate bn-1
            bn1 = calculate_bn1_log(n);
            % Calculate the non-central parameter delta
            delta = 3*sqrt(n)*C;
            
            % Calculate critical values for each alpha
            for j = 1:length(alpha_values)
                try
                    % Calculate critical value of the non-central t-distribution
                    t_crit = nctinv(1-alpha_values(j), n-1, delta);
                    % Calculate c0
                    c0 = bn1*t_crit/(3*sqrt(n));
                    % Store the result
                    results(i,j) = c0;
                catch
                    results(i,j) = NaN;
                end
            end
        catch
            results(i,:) = NaN;
        end
    end

    % Create a table
    T = array2table([n_values', results], ...
        'VariableNames', {'n', 'alpha_01', 'alpha_025', 'alpha_05'});

    % Display results
    fprintf('Critical values c0 for C = %.2f, n = 10(5)505 and α = 0.01, 0.025, 0.05\n\n', C);

    % Use fixed-width format for aligned output
    fprintf('%-8s%14s%14s%14s\n', 'n', 'α = 0.01', 'α = 0.025', 'α = 0.05');
    fprintf('%-8s%14s%14s%14s\n', '--', '--------', '--------', '--------');

    % Format and output each row
    for i = 1:height(T)
        if any(isnan(T{i,2:end}))
            % If there are NaN values, display in a special format
            fprintf('%-8d%14s%14s%14s\n', T.n(i), 'NaN', 'NaN', 'NaN');
        else
            % Display normal values with three decimal places and fixed width
            fprintf('%-8d%14.3f%14.3f%14.3f\n', ...
                T.n(i), T.alpha_01(i), T.alpha_025(i), T.alpha_05(i));
        end
    end
end

% Execute the main function
calculate_critical_values();
