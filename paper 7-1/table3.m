% 固定參數
alpha = 0.05;  % 生產者風險 α
beta = 0.05;   % 消費者風險 β
parameter_sets = [
    1.33, 1.00;
    1.50, 1.33
];
xi = 1;        % 製程偏移量 ξ = 1
m_values = 1:10; % 重複抽樣次數從1到10

% 初始化結果矩陣
results = [];

% 遍歷參數組合
for k = 1:size(parameter_sets, 1)
    CAQL = parameter_sets(k, 1);
    CLTPD = parameter_sets(k, 2);
    
    % 計算規格界限係數
    bA = 3 * CAQL + xi;
    bL = 3 * CLTPD + xi;
    
    for m = m_values
        % 定義積分函數
        S1 = @(n, c0) (1 - integral(@(t) chi2cdf((n-1)*(bA*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
            (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, bA*sqrt(n)))^m - alpha;
        S2 = @(n, c0) (1 - integral(@(t) chi2cdf((n-1)*(bL*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
            (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, bL*sqrt(n)))^m - (1-beta);
        
        % 建立聯立方程組
        system_of_equations = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
        
        % 設定初始值和求解選項
        if CAQL == 1.33
            initial_guess = [80/(m^0.5), 1.1 + 0.05*m];
        else
            initial_guess = [400/(m^0.5), 1.4 + 0.05*m];
        end
        
        options = optimoptions('fsolve', ...
            'Display', 'off', ...
            'FunctionTolerance', 1e-8, ...
            'StepTolerance', 1e-8, ...
            'MaxIterations', 1000);
        
        % 求解方程組
        solution = fsolve(system_of_equations, initial_guess, options);
        n = ceil(solution(1)); % 樣本數取整數
        c0 = solution(2);      % 臨界值
        
        % 儲存結果
        results = [results; CAQL, CLTPD, m, n, c0];
        
    end
end

% 結果轉表格顯示
result_table = array2table(results, ...
    'VariableNames', {'CAQL', 'CLTPD', 'm', 'SampleSize_n', 'CriticalValue_c0'});
disp(result_table);