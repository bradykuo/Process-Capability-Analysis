% Set parameters 
n_values = 5:5:200;  % From 5 to 200 in steps of 5
cpu_values = 0.7:0.1:3.0;  % From 0.7 to 3.0 in steps of 0.1
gamma = 0.95;

% Initialize results matrix
results = zeros(length(n_values), length(cpu_values));

% Calculate each cell in the table
for i = 1:length(n_values)
    n = n_values(i);
    for j = 1:length(cpu_values)
        cpu = cpu_values(j);
        
        % Calculate standard deviation and mean for the desired CPU
        sigma = 1;  % Assume unit standard deviation
        mu = 0;    % Assume zero mean
        USL = mu + 3*sigma*cpu;  % Calculate USL based on CPU
        
        % Calculate bn-1 using logarithms
        log_bn1 = 0.5 * log(2/(n-1)) + gammaln((n-1)/2) - gammaln((n-2)/2);
        bn1 = exp(log_bn1);
        
        % Calculate R using logarithms
        log_R = 0.5 * log(2/(n-1)) + gammaln(n/2) - gammaln((n-1)/2);
        R = exp(log_R);
        
        % Initial calculations
        t = 3 * cpu * sqrt(n)/bn1;
        z = norminv(gamma, 0, 1);
        
        % Calculate polynomial coefficients
        e = R^2 - z^2 * (1 - R^2);
        f = -2 * R * (t * R^2 - z^2 * (1 - R^2) * t);
        g = (t * R^2 - z^2 * (1 - R^2) * t)^2 - z^2 * R^2 + z^4 * (1 - R^2);
        
        % Find roots and initial delta
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
                results(i,j) = delta(k)/(3 * sqrt(n));
                break
            end
        end
    end
    % Display progress
    fprintf('Completed calculations for n = %d\n', n_values(i));
end

% Display results
fprintf('\nTable 3: Lower confidence bounds CU for CPU, gamma = 0.95\n\n');
fprintf('n\\CPU ');
fprintf('%6.1f ', cpu_values);
fprintf('\n');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');

for i = 1:length(n_values)
    fprintf('%3d ', n_values(i));
    fprintf('%6.3f ', results(i,:));
    fprintf('\n');
end
    