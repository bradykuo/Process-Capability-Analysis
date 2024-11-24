% 計算方程式(9)和(10)中積分的函數
function result = calculate_integral(n, c0, C_value, xi)
    b = 3 * C_value + abs(xi);
    % 定義被積函數
    function y = integrand(t)
        % G 是自由度為 n-1 的卡方分配累積機率函數
        G_term = chi2cdf((n - 1) * (b * sqrt(n) - t).^2 ./ (9 * n * c0^2), n-1);
        % phi 是標準常態分配的機率密度函數
        phi_term = normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n));
        y = G_term .* phi_term;
    end
    % 數值積分
    result = integral(@integrand, 0, b*sqrt(n));
end

% 針對給定的 n 和參數找出 c0 的函數
function [c0_opt, fval] = find_c0(n, C_AQL, C_LTPD, alpha, beta, xi)
    % 目標函數
    function err = objective(c0)
        eq1 = calculate_integral(n, c0, C_AQL, xi) - (1 - alpha);
        eq2 = calculate_integral(n, c0, C_LTPD, xi) - beta;
        err = eq1^2 + eq2^2;
    end
    % 尋找最佳的 c0
    [c0_opt, fval] = fminbnd(@objective, 0.8, 2.0);
end

% 給定參數
xi_values = 0:0.1:2;
CAQL_values = [1.33, 1.50, 1.67, 2.00];
CLTPD = 1.00;
alpha = 0.05;
beta = 0.05;

% 初始化陣列
n_results = zeros(length(xi_values), length(CAQL_values));
c0_results = zeros(length(xi_values), length(CAQL_values));

% 計算每個 CAQL 和 xi 組合的值
for i = 1:length(CAQL_values)
    CAQL = CAQL_values(i);
    for j = 1:length(xi_values)
        xi = xi_values(j);
        % 找出滿足兩個方程式的最小 n
        n_min = 10; % 從小的 n 開始
        while true
            [c0, fval] = find_c0(n_min, CAQL, CLTPD, alpha, beta, xi);
            % 檢查方程式是否滿足條件
            eq1 = calculate_integral(n_min, c0, CAQL, xi) >= (1 - alpha);
            eq2 = calculate_integral(n_min, c0, CLTPD, xi) <= beta;
            if eq1 && eq2
                break
            end
            n_min = n_min + 1;
            if n_min > 100
                break 
            end
        end
        % 將結果存進陣列
        n_results(j,i) = n_min;
        c0_results(j,i) = c0;
    end
end

% 建立圖形
figure('Position', [100 100 1200 500]);

% 所需樣本數 n 對 ξ 的關係圖​
subplot(1,2,1);
plot(xi_values, n_results, 'LineWidth', 2);
ylim([0 100]);
grid on;
title('Required Sample Size vs \xi');
xlabel('\xi');
ylabel('Required sample size n');
legend(cellstr(num2str(CAQL_values', 'C_{AQL}=%.2f')), 'Location', 'best');

% 臨界接受值 c₀ 對 ξ 的關係圖
subplot(1,2,2);
plot(xi_values, c0_results, 'LineWidth', 2);
ylim([1.0 1.6]);
grid on;
title('Critical Acceptance Value vs \xi');
xlabel('\xi');
ylabel('Critical value C_o');
legend(cellstr(num2str(CAQL_values', 'C_{AQL}=%.2f')), 'Location', 'best');

set(gcf, 'Color', 'white');
set(findall(gcf,'-property','FontSize'),'FontSize', 12);