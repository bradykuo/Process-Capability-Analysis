function min_n()
 
global alpha beta CAQL CLTPD xi  % 全域變數
 alpha = 0.05;     % 生產者風險
 beta = 0.05;      % 消費者風險
 CAQL = 2.00;      % 可接受品質水準
 CLTPD = 1.67;     % 可容忍品質水準
 xi = 1;           % 製程偏移量

% 設定數值求解參數
 options = optimoptions('fsolve', ...
'Display', 'iter', ...
'FunctionTolerance', 1e-8, ...
'StepTolerance', 1e-8, ...
'MaxIterations', 1000);

% 設定初始值 [n; c0]
 x0 = [10; 1];

% 執行求解
try
 [solution, fval, exitflag] = fsolve(@equations, x0, options);
if exitflag > 0
 n = ceil(solution(1));    % 樣本數取整數
 c0 = solution(2);
 
 % 驗證結果
 P_AQL = calculate_P(n, c0, CAQL, xi);
 P_LTPD = calculate_P(n, c0, CLTPD, xi);
 fprintf('最小的 n = %d\n', n);
 fprintf('對應的 c0 = %.4f\n', c0);
 fprintf('P(AQL) = %.6f (應 ≥ %.6f)\n', P_AQL, 1-alpha);
 fprintf('P(LTPD) = %.6f (應 ≤ %.6f)\n', P_LTPD, beta);

 % 檢查約束條件
if P_AQL >= 1-alpha && P_LTPD <= beta
 fprintf('找到的解滿足所有約束條件\n');
else
 fprintf('警告：解不完全滿足約束條件\n');
end
else
 fprintf('求解失敗\n');
end
catch ME
 fprintf('計算過程發生錯誤: %s\n', ME.message);
end
end

% 計算接受機率的函數
function P = calculate_P(n, c0, C_value, xi)
 b = 3 * C_value + xi;
 P = integral(@(t) chi2cdf((n - 1) * (b * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
 .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), ...
 0, b * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9);
end

% 求解方程組的函數
function F = equations(x)
global alpha beta CAQL CLTPD xi
 n = x(1);
 c0 = x(2);
 P_AQL = calculate_P(n, c0, CAQL, xi);
 P_LTPD = calculate_P(n, c0, CLTPD, xi);
 F = [P_AQL - (1 - alpha);    % 生產者風險方程式
      P_LTPD - beta];         % 消費者風險方程式
end