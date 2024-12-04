function calculatePanels()
    % 開始計時
    tic;
    
    % 定義參數範圍
    n_values = 10:5:200;              % 樣本量範圍：10 到 200，間隔 5
    c_pk_hat_values_A = 0.7:0.1:1.8;  % Cpk 估計值範圍 A ：0.7 到 1.8，間隔 0.1
    c_pk_hat_values_B = 1.9:0.1:3.0;  % Cpk 估計值範圍 B：1.9 到 3.0，間隔 0.1
    
    % 計算 Panel A 的結果
    fprintf('Calculating Panel A...\n');
    panel_A = zeros(length(n_values), length(c_pk_hat_values_A));
    for i = 1:length(c_pk_hat_values_A)
        for j = 1:length(n_values)
            panel_A(j,i) = round(calculate_lcb(n_values(j), c_pk_hat_values_A(i)), 3);
        end
    end
    
    % 計算 Panel B 的結果
    fprintf('Calculating Panel B...\n');
    panel_B = zeros(length(n_values), length(c_pk_hat_values_B));
    for i = 1:length(c_pk_hat_values_B)
        for j = 1:length(n_values)
            panel_B(j,i) = round(calculate_lcb(n_values(j), c_pk_hat_values_B(i)), 3);
        end
    end
    
    % 創建結果表格
    result_table_A = array2table([n_values' panel_A], 'VariableNames', ...
        ['n', cellstr(num2str(c_pk_hat_values_A', '%.1f'))']);
    result_table_B = array2table([n_values' panel_B], 'VariableNames', ...
        ['n', cellstr(num2str(c_pk_hat_values_B', '%.1f'))']);
    
    % 印出結果
    disp('Panel A:');
    disp(result_table_A);
    disp('Panel B:');
    disp(result_table_B);
    
    % 印出執行時間
    execution_time = toc;
    fprintf('\nTotal execution time: %.2f seconds\n', execution_time);
end

% 計算置信下界的函數
function lcb = calculate_lcb(n, Cpk_hat, xi_hat, gamma)
    % 設定默認參數值
    if nargin < 3
        xi_hat = 1.0;  % 默認 xi_hat 值
    end
    if nargin < 4
        gamma = 0.95;  % 默認置信水平
    end
    
    % 設定 fzero 函數的選項（關閉輸出顯示）
    options = optimset('Display', 'off');
    
    % 使用 fzero 尋找根
    lcb = fzero(@(C) integrate_func(C, n, Cpk_hat, xi_hat, gamma), [0, Cpk_hat], options);
end

% 積分函數
function result = integrate_func(C, n, Cpk_hat, xi_hat, gamma)
    b = 3*C + abs(xi_hat);
    
    % 定義被積函數
    integrand = @(t) G_func(t, n, b, Cpk_hat) .* phi_func(t, xi_hat, n);
    
    % 進行數值積分
    result = integral(integrand, 0, b*sqrt(n)) - (1 - gamma);
end

% 計算卡方分布的累積分布函數
function G = G_func(t, n, b, Cpk_hat)
    G = chi2cdf((n-1)*(b*sqrt(n)-t).^2 ./ (9*n*Cpk_hat^2), n-1);
end

% 計算正態分布的密度函數
function phi = phi_func(t, xi_hat, n)
    phi = normpdf(t + xi_hat*sqrt(n)) + normpdf(t - xi_hat*sqrt(n));
end

% 執行主函數
calculatePanels()