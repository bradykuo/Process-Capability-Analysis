% Set parameters
n_values = 10:10:50;           % n from 10 to 50 in steps of 10
cpu_values = 0.7:0.1:1.5;      % CPU from 0.7 to 1.5 in steps of 0.1
gamma = 0.95;                  % 95% confidence level

% Initialize results matrices
results_new = zeros(length(n_values), length(cpu_values));
results_chou = zeros(length(n_values), length(cpu_values));

% Calculate for each combination
for i = 1:length(n_values)
    n = n_values(i);
    for j = 1:length(cpu_values)
        cpu = cpu_values(j);
        
        % This paper's method (New)
        log_bn1 = 0.5 * log(2/(n-1)) + gammaln((n-1)/2) - gammaln((n-2)/2);
        bn1 = exp(log_bn1);
        log_R = 0.5 * log(2/(n-1)) + gammaln(n/2) - gammaln((n-1)/2);
        R = exp(log_R);
        
        t = 3 * cpu * sqrt(n)/bn1;
        z = norminv(gamma, 0, 1);
        
        e = R^2 - z^2 * (1 - R^2);
        f = -2 * R * (t * R^2 - z^2 * (1 - R^2) * t);
        g = (t * R^2 - z^2 * (1 - R^2) * t)^2 - z^2 * R^2 + z^4 * (1 - R^2);
        
        rt = roots([e f g]);
        x = nctinv(gamma, n-1, rt(2));
        x1 = nctinv(gamma, n-1, rt(2) + 0.01);
        x2 = (x1 - x)/9.5;
        y = abs(x - sqrt(n) * 3 * cpu/bn1);
        y1 = y/x2;
        delta1 = rt(2) + y1 * 0.001;
        
        % Iterate to find CU
        delta = delta1 + 0.001 : 0.001 : delta1 + 0.5;
        for k = 1:500
            x3 = nctinv(gamma, n-1, delta(k));
            if(abs(x3 - sqrt(n) * 3 * cpu/bn1)) <= 0.01
                results_new(i,j) = delta(k)/(3 * sqrt(n));
                break
            end
        end
        
        % Chou's method (1990)
        % Use binary search to find cu that satisfies:
        % Pr[T_{n-1}(δ = 3√n c_U) ≤ 3*CPU_hat*√n] = γ
        cu_low = 0;
        cu_high = cpu;
        tolerance = 1e-6;
        max_iter = 100;
        iter = 0;
        
        while (cu_high - cu_low) > tolerance && iter < max_iter
            cu = (cu_low + cu_high) / 2;
            
            % Calculate noncentral parameter δ = 3√n c_U
            delta = 3 * sqrt(n) * cu;
            
            % Calculate critical value = 3*CPU_hat*√n
            t_crit = 3 * cpu * sqrt(n);
            
            % Calculate probability
            prob = nctcdf(t_crit, n-1, delta);
            
            % Update bounds based on probability
            if prob < gamma
                cu_high = cu;
            else
                cu_low = cu;
            end
            
            iter = iter + 1;
        end
        
        results_chou(i,j) = (cu_low + cu_high) / 2;
    end
end

% Display results
fprintf('\nTable: Comparisons of the LCB CU between Chous and our approach (New)\n');
fprintf('CPU = 0.7(0.1)1.5, n = 10(10)50, gamma = 0.95\n\n');
fprintf('n      Method    ');
fprintf('%6.1f ', cpu_values);
fprintf('\n');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:length(n_values)
    % Print New results
    fprintf('%-4d   New     ', n_values(i));
    fprintf('%6.2f ', results_new(i,:));
    fprintf('\n');
    
    % Print Chou results
    fprintf('       Chou    ');
    fprintf('%6.2f ', results_chou(i,:));
    fprintf('\n\n');
end
