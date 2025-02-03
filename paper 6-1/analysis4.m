function analysis4()
    % Fixed parameters
    n = 50;
    ka = 1.17; % acceptance critical value
    kr1 = 1.0; % first rejection critical value
    kr2 = 0.8; % second rejection critical value
    xi = 1.0;

    % Generate Cpk values
    cpk_values = 0.5:0.01:2.0;

    % Initialize arrays
    p_single = zeros(size(cpk_values));
    p_rgs1 = zeros(size(cpk_values));
    p_rgs2 = zeros(size(cpk_values));

    % Calculate probabilities
    for i = 1:length(cpk_values)
        cpk = cpk_values(i);
        b = 3 * cpk + xi;

        Pa = compute_F_cpk(ka, n, b, xi);
        Pr1 = 1 - compute_F_cpk(kr1, n, b, xi);
        Pr2 = 1 - compute_F_cpk(kr2, n, b, xi);

        p_single(i) = Pa;
        p_rgs1(i) = Pa / (Pa + Pr1);
        p_rgs2(i) = Pa / (Pa + Pr2);
    end

    % Create plot
    figure('Position', [100 100 800 600]);
    hold on;

    plot(cpk_values, p_single, '--', 'Color', [0 0 1], 'LineWidth', 1.5);
    plot(cpk_values, p_rgs1, '-', 'Color', [0 0.5 0], 'LineWidth', 1.5);
    plot(cpk_values, p_rgs2, '-', 'Color', [1 0 0], 'LineWidth', 1.5);

    grid on;
    xlabel('Cpk Value');
    ylabel('Probability of Acceptance');
    title('OC Curves of Sampling Plans n = 50');
    legend('Variables Single Sampling Plan', ...
           'Variables RGS Plan (kr = 1)', ...
           'Variables RGS Plan (kr = 0.8)', ...
           'Location', 'southeast');

    xlim([0.5 2.0]);
    ylim([0 1.0]);
    grid minor;
    hold off;
end

function F = compute_F_cpk(ka, n, b, xi)
    function y = integrand(t)
        chi_term = (n-1) * (b*sqrt(n) - t).^2 / (9*n*ka^2);
        phi_terms = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = chi2cdf(chi_term, n-1) .* phi_terms;
    end
    F = integral(@integrand, 0, b*sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9);
end