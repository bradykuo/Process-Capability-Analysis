% Define parameters 
xi_values = 0:0.05:3.0;  % ξ values from 0 to 3.0 
beta = 0.05; 
alpha = 0.05; 
CLTPD = 1.00; 

% Different CAQL values 
CAQL_values = [1.33, 1.50, 1.67, 2.00]; 

% Initialize arrays to store sample sizes 
n_values = zeros(length(CAQL_values), length(xi_values)); 

% Calculate sample sizes for each CAQL and ξ value 
for i = 1:length(CAQL_values)
    CAQL = CAQL_values(i);
    for j = 1:length(xi_values)
        xi = xi_values(j);
        
        % Calculate b1 and b2 
        b1 = 3 * CAQL * (1 + xi^2)^0.5;
        b2 = 3 * CLTPD * (1 + xi^2)^0.5;
        
        % Adjust initial guesses based on CAQL value
        if CAQL >= 1.67
            if j == 1  
                n = 20;  % Start with smaller n for higher CAQL
                C0 = 1.4;  % Higher initial C0 for higher CAQL
            else
                % Use previous successful solution as initial guess
                if ~isnan(n_values(i,j-1))
                    n = n_values(i,j-1);
                    C0 = 1.4;
                else
                    n = 20;
                    C0 = 1.4;
                end
            end
        else
            n = 50;
            C0 = 1.2;
        end
        
        % Define S1 and S2 according to equations (10) and (11) 
        S1 = @(n, C0) integral(@(t) chi2cdf((b1.^2.*n)./(9.*C0.^2) - t.^2, n-1) .* ...
            (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
            0, (b1.*sqrt(n))./(3.*C0)) - (1-alpha);
        
        S2 = @(n, C0) integral(@(t) chi2cdf((b2.^2.*n)./(9.*C0.^2) - t.^2, n-1) .* ...
            (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
            0, (b2.*sqrt(n))./(3.*C0)) - beta;
        
        % Define objective function
        obj_fun = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
        
        % Try multiple initial guesses
        initial_guesses = [
            n, C0;
            n*1.2, C0*1.1;
            n*0.8, C0*0.9;
        ];
        
        options = optimset('Display', 'off');
        success = false;
        
        for guess = 1:size(initial_guesses, 1)
            try
                x0 = initial_guesses(guess,:);
                x_sol = fsolve(obj_fun, x0, options);
                
                % Check if solution is reasonable
                if x_sol(1) > 0 && x_sol(1) < 200 && x_sol(2) > 0.8 && x_sol(2) < 2.0
                    n_values(i,j) = x_sol(1);
                    success = true;
                    break;
                end
            catch
                continue;
            end
        end
        
        if ~success
            n_values(i,j) = NaN;
        end
    end
end

% Create the plot 
figure;
hold on;

% Plot lines for each CAQL value 
for i = 1:length(CAQL_values)
    % Filter out NaN values and interpolate
    valid_idx = ~isnan(n_values(i,:));
    if any(valid_idx)
        plot(xi_values(valid_idx), n_values(i,valid_idx), 'LineWidth', 1.5);
    end
end

% Customize the plot 
xlabel('\xi');
ylabel('sample size n');
axis([0 2.0 0 80]);
box on;

% Add legend 
legend_entries = cell(length(CAQL_values), 1);
for i = 1:length(CAQL_values)
    legend_entries{i} = sprintf('C_{AQL} = %.2f', CAQL_values(i));
end
legend(legend_entries, 'Location', 'northeast');
hold off;