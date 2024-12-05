% 計算抽樣計畫的函數
function [n_results, c0_results] = calculate_sampling_plan(xi_values, CAQL_values, CLTPD, alpha, beta)

% 定義積分計算函數
function [result] = calculate_integral(n, c0, C_value, xi)
 b = 3 * C_value + abs(xi);
function y = integrand(t)
 % 計算卡方分配和常態分配項
 G_term = chi2cdf((n - 1) * (b * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
 phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
 y = G_term .* phi_term;
end
 result = integral(@integrand, 0, b*sqrt(n));
end

% 定義 S1 和 S2 聯立方程式
function F = equations(x, CAQL, CLTPD, alpha, beta, xi)
 n = x(1);
 c0 = x(2);
 S1 = calculate_integral(n, c0, CAQL, xi) - (1 - alpha);
 S2 = calculate_integral(n, c0, CLTPD, xi) - beta;
 F = [S1; S2];
end

n_results = zeros(length(xi_values), length(CAQL_values));
c0_results = zeros(length(xi_values), length(CAQL_values));

% 對每個參數組合求解
for i = 1:length(CAQL_values)
 CAQL = CAQL_values(i);
for j = 1:length(xi_values)
 xi = xi_values(j);
 
 % 設定求解參數
 obj_fun = @(x) equations(x, CAQL, CLTPD, alpha, beta, xi);
 x0 = [50, 1.2]; % 初始猜測值
 options = optimoptions('fsolve', 'Display', 'off');
 
 % 使用 fsolve 求解方程組
try
 [x_sol, fval, exitflag] = fsolve(obj_fun, x0, options);
if exitflag > 0
 n_results(j,i) = ceil(x_sol(1)); % 樣本數取整數
 c0_results(j,i) = x_sol(2);
else
 n_results(j,i) = NaN;
 c0_results(j,i) = NaN;
end
catch
 n_results(j,i) = NaN;
 c0_results(j,i) = NaN;
end
end
end
end

xi_values = 0:0.1:2;         % 製程偏移量
CAQL_values = [1.33, 1.50, 1.67, 2.00];  % 可接受品質水準
CLTPD = 1.00;                % 可容忍品質水準
alpha = 0.05;                % 生產者風險
beta = 0.05;                 % 消費者風險

% 執行計算
[n_results, c0_results] = calculate_sampling_plan(xi_values, CAQL_values, CLTPD, alpha, beta);

% 繪圖設定
figure('Position', [100 100 1200 500]);

% 繪製樣本數與 ξ 關係圖
subplot(1,2,1);
plot(xi_values, n_results, 'LineWidth', 2);
ylim([0 100]);
grid on;
title('Required Sample Size vs \xi');
xlabel('\xi');
ylabel('Required sample size n');
legend(cellstr(num2str(CAQL_values', 'C_{AQL}=%.2f')), 'Location', 'best');

% 繪製臨界值與 ξ 關係圖
subplot(1,2,2);
plot(xi_values, c0_results, 'LineWidth', 2);
ylim([1.0 1.6]);
grid on;
title('Critical Acceptance Value vs \xi');
xlabel('\xi');
ylabel('Critical value C_o');
legend(cellstr(num2str(CAQL_values', 'C_{AQL}=%.2f')), 'Location', 'best');

% 設定圖形外觀
set(gcf, 'Color', 'white');
set(findall(gcf,'-property','FontSize'),'FontSize', 12);