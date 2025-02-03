% Function to calculate bn-1
function b = calculate_bn1(n)
    b = sqrt(2/(n-1)) * gamma((n-1)/2) / gamma((n-2)/2);
end

% Function to calculate the integrand
function val = integrand(y, x, n, delta, bn1)
    term1 = y.^((n-2)/2);
    term2 = exp(-0.5 * (y + (3*x*sqrt(n*y)/(bn1*sqrt(n-1)) - delta).^2));
    val = term1 .* term2;
end

% Set parameters
CPU_values = [1.00, 1.25, 1.45, 1.60];
sample_sizes = [10, 20, 30, 40, 50];
x = linspace(0, 3, 300); 

% Create figure with subplots
figure('Position', [100, 100, 1200, 800]);
subplot_titles = {'(a) CPU = 1.00', '(b) CPU = 1.25', '(c) CPU = 1.45', '(d) CPU = 1.60'};

for subplot_idx = 1:length(CPU_values)
    subplot(2, 2, subplot_idx);
    hold on;
    
    CPU = CPU_values(subplot_idx);
    
    % Plot for each sample size
    for n = sample_sizes
        % Calculate constants
        bn1 = calculate_bn1(n);
        delta = 3*sqrt(n)*CPU;
        
        % Calculate PDF values
        pdf_values = zeros(size(x));
        for i = 1:length(x)
            % Numerical integration
            integral_val = integral(@(y) integrand(y, x(i), n, delta, bn1), 0, Inf);
            
            % Calculate coefficient
            coef = (3*sqrt(n/(n-1))*2^(-n/2))/(bn1*sqrt(pi)*gamma((n-1)/2));
            
            pdf_values(i) = coef * integral_val;
        end
        
        % Plot the PDF
        plot(x, pdf_values, 'LineWidth', 1.5);
    end
    
    % Customize subplot
    grid on;
    xlabel('x');
    ylabel('f(x)');
    title(subplot_titles{subplot_idx});
    ylim([0 3.5]);
    xlim([0 3]);
    
    % Add legend only to the first subplot
    if subplot_idx == 1
        legend('n=10', 'n=20', 'n=30', 'n=40', 'n=50', 'Location', 'northeast');
    end
end

% Adjust spacing between subplots
set(gcf, 'Color', 'white');
sgtitle('PDF Plots of C_{PU} for Different CPU Values', 'FontSize', 14);

% Make the figure more compact
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.8, 0.8]);