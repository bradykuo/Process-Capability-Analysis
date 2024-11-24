% 初始化參數
alpha_values = [0.01, 0.025, 0.05, 0.075, 0.10]; % 供應商風險 α
beta_values = [0.01, 0.025, 0.05, 0.075, 0.10]; % 客戶風險 β
parameter_sets = [
    1.33, 1.00;  % CAQL, CLTPD pair 1
    1.50, 1.33;  % CAQL, CLTPD pair 2
    1.67, 1.33;  % CAQL, CLTPD pair 3
    2.00, 1.67   % CAQL, CLTPD pair 4
];
xi = 1; % 設定 ξ = 1

% 儲存結果
results = [];

% 遍歷所有組合
for k = 1:size(parameter_sets, 1)
    CAQL = parameter_sets(k, 1); %k組的第一個值
    CLTPD = parameter_sets(k, 2);  %k組的第二個值  
    for i = 1:length(alpha_values)
        alpha = alpha_values(i);
        for j = 1:length(beta_values)
            beta = beta_values(j);
            
            % 計算 b1 和 b2
            b1 = 3 * CAQL + xi;
            b2 = 3 * CLTPD + xi;

            % 定義 S1 和 S2 函數的積分
            S1 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
                .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b1 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1 - alpha);

            S2 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
                .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b2 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;

            % 定義方程系統，求解聯立方程組S1=0、S2=0
            system_of_equations = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];

            % 設定初始猜測值和求解選項
            initial_guess = [50 + 20 * (CAQL - 1), 1.0 + 0.5 * (CAQL - 1)];%根據CAQL調整初始猜測值
            options = optimoptions('fsolve', ...%用fsolve求解
                'Display', 'off', ...
                'FunctionTolerance', 1e-8, ...
                'StepTolerance', 1e-8, ...
                'MaxIterations', 1000);

            % 求解
            try
                solution = fsolve(system_of_equations, initial_guess, options);
                n = ceil(solution(1));  %向上取整數
                c0 = solution(2);
            catch
                n = NaN;
                c0 = NaN;
            end
            % 儲存結果
            results = [results; alpha, beta, CAQL, CLTPD, n, c0];
            
            % 顯示進度
            fprintf('Processing α=%.3f, β=%.3f, CAQL=%.2f, CLTPD=%.2f: n=%.0f, c0=%.4f\n', ...
                alpha, beta, CAQL, CLTPD, n, c0);
        end
    end
end

% 將結果轉換為表格
result_table = array2table(results, ...
    'VariableNames', {'Alpha', 'Beta', 'CAQL', 'CLTPD', 'SampleSize_n', 'CriticalValue_c0'});
disp(result_table);
