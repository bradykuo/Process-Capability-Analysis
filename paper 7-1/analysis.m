alpha_values = 0.05; % 生產者風險 α
beta_values = 0.05; % 消費者風險 β
parameter_sets = [
 1.33, 1.00; 
 1.50, 1.00; 
 1.67, 1.00;
 2.00, 1.00
];
xi = 1;  % 製程偏移量 ξ = 1
m = 2;   % 重複抽樣次數

% 初始化結果矩陣
results = [];

% 遍歷參數組合
for k = 1:size(parameter_sets, 1)
 CAQL = parameter_sets(k, 1);
 CLTPD = parameter_sets(k, 2);
    for i = 1:length(alpha_values)
     alpha = alpha_values(i);
        for j = 1:length(beta_values)
         beta = beta_values(j);
        
         % 計算規格界限係數
         bA = 3 * CAQL + xi;
         bL = 3 * CLTPD + xi;
        
         % 定義積分函數
         S1 = @(n, c0) (1 - integral(@(t) chi2cdf((n-1)*(bA*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
            (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, bA*sqrt(n)))^m - (alpha);

         S2 = @(n, c0) (1 - integral(@(t) chi2cdf((n-1)*(bL*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
            (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, bL*sqrt(n)))^m - (1-beta);

         % 建立聯立方程組
         system_of_equations = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
        
         % 設定初始值和求解選項
         initial_guess = [10, 1];
         options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxIterations', 1000);
        
         % 求解方程組
         solution = fsolve(system_of_equations, initial_guess, options);
         n = ceil(solution(1));  % 樣本數取整數
         c0 = solution(2);       % 臨界值
         
         % 儲存結果
         results = [results; alpha, beta, CAQL, CLTPD, n, c0];
        
        end
    end
end

% 結果轉表格顯示
result_table = array2table(results, ...
'VariableNames', {'Alpha', 'Beta', 'CAQL', 'CLTPD', 'SampleSize_n', 'CriticalValue_c0'});
disp(result_table);