% 定義計算 c0 的主函數
function c0 = calculate_c0(C, Cp, n, alpha)
    % 定義被積分函數
    function result = integrand(y, c0)
        % 計算卡方分布累積分布函數項
        G_term = chi2cdf((n - 1) * (3 * Cp * sqrt(n) - y).^2 ./ (9 * n * c0^2), n - 1);
        % 計算正態分布密度函數項
        f_term = normpdf(y + 3 * (Cp - C) * sqrt(n)) + normpdf(y - 3 * (Cp - C) * sqrt(n));

        result = G_term .* f_term;
    end

    % 定義目標函數（用於尋找根）
    function result = objective(c0)
        % 計算定積分
        integral_result = integral(@(y) integrand(y, c0), 0, 3 * Cp * sqrt(n));
        % 計算與目標值 alpha 的差
        result = integral_result - alpha;
    end
    
    % 設定初始值
    initial_value = 1;
    
    % 尋找適當的區間進行根查找
    % lower = 0.5;  % 下界
    % upper = 5;    % 上界
    % 當區間端點函數值乘積大於 0 時，調整區間
    % while objective(lower) * objective(upper) > 0
    %     if objective(lower) > 0
    %         lower = lower / 2;  % 縮小下界
    %     else
    %         upper = upper * 2;  % 擴大上界
    %     end
    % end
    % [lower, upper]
    
    % 使用 fzero 函數進行根查找
    options = optimset('TolX', 1e-8);  % 設定收斂容差
    c0 = fzero(@objective, initial_value, options);  % 尋找根
end

% 設定參數值
C_values = [1.00, 1.33, 1.5, 1.67, 2.00];  % C 值
n_values = 10:5:405;                        % 樣本數序列（從 10 到 405，間隔 5 ）
alpha_values = [0.01, 0.025, 0.05];         % 顯著性水準

% 創建儲存結果的陣列
num_n = length(n_values);          
num_alpha = length(alpha_values);  
num_C = length(C_values);          

% 對每個 C 值進行循環計算
for c_idx = 1:num_C
    C = C_values(c_idx);
    
    % 創建儲存當前 C 值結果的矩陣
    results = zeros(num_n, num_alpha);
    
    for i = 1:num_n
        n = n_values(i);
        % 根據樣本數決定 Cp 值
        if n < 100
            Cp = C + 0.33;  % 樣本數小於 100 時的調整
        else
            Cp = C + 0.12;  % 樣本數大於等於 100 時的調整
        end
        
        % 計算不同顯著性水準的結果
        for j = 1:num_alpha
            alpha = alpha_values(j);
            results(i, j) = calculate_c0(C, Cp, n, alpha);
        end
    end
    
    % 輸出結果表格
    fprintf('\nCritical values c0 for C = %.2f, n = 10(5)230 and α = 0.01, 0.025, 0.05\n', C);
    fprintf('---------------------------------------------------------------\n');
    fprintf('     n    α = 0.01   α = 0.025    α = 0.05\n');
    fprintf('---------------------------------------------------------------\n');
    
    % 列印結果數據
    for i = 1:num_n
        fprintf('%6d', n_values(i));  
        for j = 1:num_alpha
            fprintf('%11.3f', results(i, j));  % 精確到小數點後3位
        end
        fprintf('\n');  
    end
    fprintf('---------------------------------------------------------------\n\n');
end