% Define parameters
alpha_values = 0.01:0.01:0.10;
beta_values = 0.01:0.01:0.10;
[BETA, ALPHA] = meshgrid(beta_values, alpha_values);  % Swapped order here
xi = 0;  % Set ξ = 0 (according to the paper)

% Initialize result matrices
N1 = zeros(size(BETA));
C01 = zeros(size(BETA));
N2 = zeros(size(BETA));
C02 = zeros(size(BETA));

% Calculate for (CAQL, CLTPD) = (1.05, 1.00)
CAQL1 = 1.05;
CLTPD1 = 1.00;
b1_1 = 3 * CAQL1 * (1 + xi^2)^0.5;
b2_1 = 3 * CLTPD1 * (1 + xi^2)^0.5;

% Calculate for (CAQL, CLTPD) = (1.10, 1.00)
CAQL2 = 1.10;
CLTPD2 = 1.00;
b1_2 = 3 * CAQL2 * (1 + xi^2)^0.5;
b2_2 = 3 * CLTPD2 * (1 + xi^2)^0.5;

% Calculate values for both sets
for i = 1:length(alpha_values)
    for j = 1:length(beta_values)
        % For (1.05, 1.00)
        S1_1 = @(n, c0) integral(@(t) chi2cdf((b1_1.^2.*n)./(9.*c0.^2) - t.^2, n-1) .* ...
            (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
            0, b1_1.*sqrt(n)./(3.*c0)) - (1-ALPHA(i,j));
        
        S2_1 = @(n, c0) integral(@(t) chi2cdf((b2_1.^2.*n)./(9.*c0.^2) - t.^2, n-1) .* ...
            (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
            0, b2_1.*sqrt(n)./(3.*c0)) - BETA(i,j);
        
        % For (1.10, 1.00)
        S1_2 = @(n, c0) integral(@(t) chi2cdf((b1_2.^2.*n)./(9.*c0.^2) - t.^2, n-1) .* ...
            (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
            0, b1_2.*sqrt(n)./(3.*c0)) - (1-ALPHA(i,j));
        
        S2_2 = @(n, c0) integral(@(t) chi2cdf((b2_2.^2.*n)./(9.*c0.^2) - t.^2, n-1) .* ...
            (normpdf(t + xi.*sqrt(n)) + normpdf(t - xi.*sqrt(n))), ...
            0, b2_2.*sqrt(n)./(3.*c0)) - BETA(i,j);
        
        % Solve equations
        try
            % For (1.05, 1.00)
            sys1 = @(x) [S1_1(x(1),x(2)); S2_1(x(1),x(2))];
            sol1 = fsolve(sys1, [50, 1.02], optimset('Display','off'));
            N1(i,j) = sol1(1);
            C01(i,j) = sol1(2);
            
            % For (1.10, 1.00)
            sys2 = @(x) [S1_2(x(1),x(2)); S2_2(x(1),x(2))];
            sol2 = fsolve(sys2, [50, 1.04], optimset('Display','off'));
            N2(i,j) = sol2(1);
            C02(i,j) = sol2(2);
        catch
            N1(i,j) = NaN;
            C01(i,j) = NaN;
            N2(i,j) = NaN;
            C02(i,j) = NaN;
        end
    end
end

% Create subplots
figure('Position', [100 100 800 1000]);

% Subplot (a)
subplot(4,1,1)
surf(BETA, ALPHA, C01, 'EdgeColor', 'none');
ylabel('\alpha');
xlabel('\beta');
zlabel('C_0');
title('(a) C_0 for (C_{AQL}, C_{LTPD}) = (1.05, 1.00)');
zlim([1.01 1.04]);
view(45, 35);

% Subplot (b)
subplot(4,1,2)
surf(BETA, ALPHA, N1, 'EdgeColor', 'none');
ylabel('\alpha');
xlabel('\beta');
zlabel('n');
title('(b) n for (C_{AQL}, C_{LTPD}) = (1.05, 1.00)');
zlim([0 4000]);
view(45, 35);

% Subplot (c)
subplot(4,1,3)
surf(BETA, ALPHA, C02, 'EdgeColor', 'none');
ylabel('\alpha');
xlabel('\beta');
zlabel('C_0');
title('(c) C_0 for (C_{AQL}, C_{LTPD}) = (1.10, 1.00)');
zlim([1.02 1.08]);
view(45, 35);

% Subplot (d)
subplot(4,1,4)
surf(BETA, ALPHA, N2, 'EdgeColor', 'none');
ylabel('\alpha');
xlabel('\beta');
zlabel('n');
title('(d) n for (C_{AQL}, C_{LTPD}) = (1.10, 1.00)');
zlim([0 1500]);
view(45, 35);

% Adjust overall appearance
colormap('jet');
set(gcf, 'Color', 'white');