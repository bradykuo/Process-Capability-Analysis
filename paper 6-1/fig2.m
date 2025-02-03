function fig2()
    % Fixed parameters
    n = 100;  % sample size
    ka = 1.214;  % acceptance critical value
    kr = 1.05;  % rejection critical value
    xi = 1.0;
    
    % Generate Cpk values
    cpk_values = 0.5:0.01:2.0;
    
    % Initialize arrays to store probabilities
    p_single = zeros(size(cpk_values));
    p_rgs = zeros(size(cpk_values));
    
    % Calculate probabilities for each Cpk value
    for i = 1:length(cpk_values)
        cpk = cpk_values(i);
        b = 3 * cpk + xi;
        
        % Calculate Pa and Pr
        Pa = compute_F_cpk(ka, n, b, xi);
        Pr = 1 - compute_F_cpk(kr, n, b, xi);
        
        % Single sampling plan
        p_single(i) = Pa;
        
        % RGS plan
        p_rgs(i) = Pa / (Pa + Pr);
    end
    
    % Create the plot
    figure('Position', [100 100 800 600]);
    hold on;
    
    % Plot single sampling plan (blue dashed line)
    plot(cpk_values, p_single, '--', 'Color', [0 0 1], 'LineWidth', 1.5);
    
    % Plot RGS plan (green solid line)
    plot(cpk_values, p_rgs, '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
    
    % Customize the plot
    grid on;
    xlabel('Cpk Value');
    ylabel('Probability of acceptance');
    title('OC curves of a variables single sampling plan and a variables RGS plan with n = 100');
    legend('Variables Single Sampling Plan', 'Variables RGS Plan', 'Location', 'southeast');
    
    % Set axis limits
    xlim([0.5 2.0]);
    ylim([0 1.0]);
    
    % Add minor grid
    grid minor;
    
    hold off;
end

function F = compute_F_cpk(ka, n, b, xi)
    % Define integrand function
    function y = integrand(t)
        chi_term = (n-1) * (b*sqrt(n) - t).^2 / (9*n*ka^2);
        phi_terms = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = chi2cdf(chi_term, n-1) .* phi_terms;
    end
    
    % Compute integral
    F = integral(@integrand, 0, b*sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9);
end