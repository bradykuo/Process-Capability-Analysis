function simulateCpkTable()
    % Parameters
    rng(123);  % for reproducibility
    Cpk_values = [1, 1.33, 1.5, 2, 2.5, 3];
    n_values = [10, 20, 30, 40, 50];
    alpha_values = [0.1, 0.05];
    num_outer_simulations = 10000;
    num_inner_simulations = 10000;
    L = 7;
    U = 14;
    mu = 10;
    d = (U - L) / 2;
    M = (U + L) / 2;

    % Pre-allocate results array
    num_Cpk = length(Cpk_values);
    num_n = length(n_values);
    num_alpha = length(alpha_values);
    results = cell(num_Cpk, 1);
    
    % Start timing
    tic;
    
    % Main simulation loop
    parfor i = 1:num_Cpk
        Cpk = Cpk_values(i);
        sigma = (d - abs(mu - M)) / (3 * Cpk);
        temp_results = cell(num_n, num_alpha);
        
        for j = 1:num_n
            n = n_values(j);
            for k = 1:num_alpha
                alpha = alpha_values(k);
                fprintf('Processing Cpk = %.2f, n = %d, alpha = %.2f\n', Cpk, n, alpha);
                sim_result = simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, ...
                    num_outer_simulations, num_inner_simulations);
                temp_results{j,k} = sim_result;
            end
        end
        results{i} = temp_results;
    end
    
    % Format and display results
    formatted_results = format_results(results, Cpk_values, n_values, alpha_values);
    
    % Display execution time
    execution_time = toc;
    fprintf('\nTotal execution time: %.2f seconds\n', execution_time);
end

function result = simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, num_outer_simulations, num_inner_simulations)
    % Generate random samples
    x = normrnd(mu, sigma, [num_outer_simulations, n]);
    x_bar = mean(x, 2);
    s = std(x, 0, 2);
    s2 = s.^2;
    
    % Calculate Gpk values
    Gpk_values = zeros(num_outer_simulations, 1);
    for i = 1:num_outer_simulations
        Gpk_values(i) = quantile(calculate_Gpk(x_bar(i), s2(i), n, L, U, num_inner_simulations), alpha);
    end
    
    % Calculate other limits
    other_limits = zeros(num_outer_simulations, 4);
    for i = 1:num_outer_simulations
        other_limits(i,:) = calculate_other_limits(x_bar(i), s(i), n, L, U, alpha);
    end
    
    % Combine limits
    limits = [Gpk_values, other_limits];
    
    % Calculate coverage and expected values
    coverage = mean(limits <= Cpk, 1);
    expected_values = mean(limits, 1);
    
    result = struct('coverage', round(coverage, 4), ...
                   'expected_values', round(expected_values, 4));
end

function Gpk = calculate_Gpk(x_bar, s2, n, L, U, num_simulations)
    d = (U - L) / 2;
    M = (U + L) / 2;
    
    Z = randn(num_simulations, 1);
    U2 = chi2rnd(n-1, [num_simulations, 1]);
    
    T_mu = x_bar - sqrt((n-1)/n) * (Z * sqrt(s2) / sqrt(n));
    T_sigma2 = s2 * (n-1) ./ U2;
    
    Gpk = (d - abs(T_mu - M)) ./ (3 * sqrt(T_sigma2));
end

function limits = calculate_other_limits(x_bar, s, n, L, U, alpha)
    Cpk_hat = calculate_Cpk_hat(x_bar, s, L, U);
    z = norminv(1 - alpha);
    
    % Calculate limits
    Bpk = Cpk_hat - z * sqrt(1/(9*n) + Cpk_hat^2/(2*(n-1)));
    KHpk = Cpk_hat * (1 - z/sqrt(2*(n-1)));
    Hpk = Cpk_hat - z * sqrt((n-1)/(9*n*(n-3)) + Cpk_hat^2*(1/(2*(n-3))*(1 + 6/(n-1))));
    Npk = sqrt(1-2/(5*(n-1)))*Cpk_hat - z*sqrt(Cpk_hat^2/(2*(n-1)) + 1/(9*n));
    
    limits = [Bpk, KHpk, Hpk, Npk];
end

function Cpk_hat = calculate_Cpk_hat(x_bar, s, L, U)
    Cpk_hat = min((U - x_bar)/(3*s), (x_bar - L)/(3*s));
end

function formatted_table = format_results(results, Cpk_values, n_values, alpha_values)
    % Create table headers
    headers = {'Cpk', 'n', '1-α', 'Gpk', 'Bpk', 'KHpk', 'Hpk', 'Npk', ...
              'E(Gpk)', 'E(Bpk)', 'E(KHpk)', 'E(Hpk)', 'E(Npk)'};
    
    % Initialize data matrix
    num_rows = length(Cpk_values) * length(n_values) * length(alpha_values);
    data = zeros(num_rows, length(headers));
    row = 1;
    
    % Fill data matrix
    for i = 1:length(Cpk_values)
        for j = 1:length(n_values)
            for k = 1:length(alpha_values)
                temp_result = results{i}{j,k};
                data(row,:) = [Cpk_values(i), n_values(j), 1-alpha_values(k), ...
                              temp_result.coverage, temp_result.expected_values];
                row = row + 1;
            end
        end
    end
    
    % Create and format table
    formatted_table = array2table(data, 'VariableNames', headers);
    
    % Display table with headers
    disp('-------------------------------------------------------------------------------------------------------------------------------');
    disp('                                     Coverage probability                           Expected value');
    disp('-------------------------------------------------------------------------------------------------------------------------------');
    disp(formatted_table);
end

simulateCpkTable()