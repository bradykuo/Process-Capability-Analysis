function analysis2()
    % Fixed parameters
    n = 50; % sample size
    ka1 = 1.17; % first acceptance critical value
    ka2 = 1.05;  % second acceptance critical value 
    kr = 1.0;   % rejection critical value
    xi = 1.0;

    % Generate Cpk values
    cpk_values = 0.5:0.01:2.0;

    % Initialize arrays to store probabilities
    p_single = zeros(size(cpk_values));
    p_rgs1 = zeros(size(cpk_values));
    p_rgs2 = zeros(size(cpk_values));

    % Calculate probabilities for each Cpk value
    for i = 1:length(cpk_values)
        cpk = cpk_values(i);
        b = 3 * cpk + xi;

        % Calculate Pa and Pr for ka=1.17
        Pa1 = compute_F_cpk(ka1, n, b, xi);
        Pr1 = 1 - compute_F_cpk(kr, n, b, xi);

        % Calculate Pa and Pr for ka=2
        Pa2 = compute_F_cpk(ka2, n, b, xi);
        Pr2 = 1 - compute_F_cpk(kr, n, b, xi);

        % Single sampling plan
        p_single(i) = Pa1;

        % RGS plans
        p_rgs1(i) = Pa1 / (Pa1 + Pr1);
        p_rgs2(i) = Pa2 / (Pa2 + Pr2);
    end

    % Create the plot
    figure('Position', [100 100 800 600]);
    hold on;

    % Plot single sampling plan (blue dashed line)
    plot(cpk_values, p_single, '--', 'Color', [0 0 1], 'LineWidth', 1.5);

    % Plot RGS plan ka=1.17 (green solid line)
    plot(cpk_values, p_rgs1, '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);

    % Plot RGS plan ka=2 (red solid line)
    plot(cpk_values, p_rgs2, '-', 'Color', [1 0 0], 'LineWidth', 1.5);

    % Customize the plot
    grid on;
    xlabel('Cpk Value');
    ylabel('Probability of Acceptance');
    title('OC Curves of Sampling Plans n = 50');
    legend('Variables Single Sampling Plan', ...
           'Variables RGS Plan (ka = 1.17)', ...
           'Variables RGS Plan (ka = 1.05)', ...
           'Location', 'southeast');

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