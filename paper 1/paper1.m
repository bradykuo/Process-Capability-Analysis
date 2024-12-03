function simulateCpkTable()
    % 設定基本參數
    rng(123);  % 設定隨機種子序確保結果可重現
    Cpk_values = [1, 1.33, 1.5, 2, 2.5, 3];  % Cpk 值的範圍
    n_values = [10, 20, 30, 40, 50];         % 樣本數的範圍
    alpha_values = [0.1, 0.05];              % 顯著水準 alpha
    num_outer_simulations = 10000;           % 外層模擬次數
    num_inner_simulations = 10000;           % 內層模擬次數
    L = 7;                                   % 規格下限
    U = 14;                                  % 規格上限
    mu = 10;                                 % 製程平均值
    d = (U - L) / 2;                         % 規格半寬
    M = (U + L) / 2;                         % 規格中點

    % 預先配置結果陣列空間
    num_Cpk = length(Cpk_values);
    num_n = length(n_values);
    num_alpha = length(alpha_values);
    results = cell(num_Cpk, 1);
    
    % 開始計時
    tic;
    
    % 主要模擬迴圈 (使用平行運算)
    parfor i = 1:num_Cpk
        Cpk = Cpk_values(i);
        sigma = (d - abs(mu - M)) / (3 * Cpk);  % 計算對應的標準差
        temp_results = cell(num_n, num_alpha);
        
        for j = 1:num_n
            n = n_values(j);
            for k = 1:num_alpha
                alpha = alpha_values(k);
                fprintf('Processing Cpk = %.2f, n = %d, alpha = %.2f\n', Cpk, n, alpha);
                sim_result = simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, ...
                    num_outer_simulations, num_inner_simulations);
                temp_results{j,k} = sim_result;
            end
        end
        results{i} = temp_results;
    end
    
    % 格式化並顯示結果
    formatted_results = format_results(results, Cpk_values, n_values, alpha_values);
    
    % 顯示執行時間
    execution_time = toc;
    fprintf('\nTotal execution time: %.2f seconds\n', execution_time);
end

function result = simulate_Cpk(n, mu, sigma, L, U, Cpk, alpha, num_outer_simulations, num_inner_simulations)
    % 生成隨機樣本
    x = normrnd(mu, sigma, [num_outer_simulations, n]);
    x_bar = mean(x, 2);       % 計算樣本平均值
    s = std(x, 0, 2);         % 計算樣本標準差
    s2 = s.^2;                % 計算樣本變異數
    
    % 計算 Gpk 值
    Gpk_values = zeros(num_outer_simulations, 1);
    for i = 1:num_outer_simulations
        Gpk_values(i) = quantile(calculate_Gpk(x_bar(i), s2(i), n, L, U, num_inner_simulations), alpha);
    end
    
    % 計算其他指標界限
    other_limits = zeros(num_outer_simulations, 4);
    for i = 1:num_outer_simulations
        other_limits(i,:) = calculate_other_limits(x_bar(i), s(i), n, L, U, alpha);
    end
    
    % 合併所有界限
    limits = [Gpk_values, other_limits];
    
    % 計算覆蓋機率和期望值
    coverage = mean(limits <= Cpk, 1);
    expected_values = mean(limits, 1);
    
    % 回傳結構化結果
    result = struct('coverage', round(coverage, 4), ...
                   'expected_values', round(expected_values, 4));
end

function Gpk = calculate_Gpk(x_bar, s2, n, L, U, num_simulations)
    % 計算規格相關參數
    d = (U - L) / 2;
    M = (U + L) / 2;
    
    % 生成隨機數
    Z = randn(num_simulations, 1);              % 標準常態分配
    U2 = chi2rnd(n-1, [num_simulations, 1]);    % 卡方分配
    
    % 計算 T 統計量
    T_mu = x_bar - sqrt((n-1)/n) * (Z * sqrt(s2) / sqrt(n));
    T_sigma2 = s2 * (n-1) ./ U2;
    
    % 計算 Gpk 值
    Gpk = (d - abs(T_mu - M)) ./ (3 * sqrt(T_sigma2));
end

function limits = calculate_other_limits(x_bar, s, n, L, U, alpha)
    % 計算 Cpk 估計值
    Cpk_hat = calculate_Cpk_hat(x_bar, s, L, U);
    z = norminv(1 - alpha);   % 計算標準常態分配的分位數
    
    % 計算各種能力指標界限
    Bpk = Cpk_hat - z * sqrt(1/(9*n) + Cpk_hat^2/(2*(n-1)));
    KHpk = Cpk_hat * (1 - z/sqrt(2*(n-1)));
    Hpk = Cpk_hat - z * sqrt((n-1)/(9*n*(n-3)) + Cpk_hat^2*(1/(2*(n-3))*(1 + 6/(n-1))));
    Npk = sqrt(1-2/(5*(n-1)))*Cpk_hat - z*sqrt(Cpk_hat^2/(2*(n-1)) + 1/(9*n));
    
    limits = [Bpk, KHpk, Hpk, Npk];
end

function Cpk_hat = calculate_Cpk_hat(x_bar, s, L, U)
    % 計算 Cpk 估計值
    Cpk_hat = min((U - x_bar)/(3*s), (x_bar - L)/(3*s));
end

function formatted_table = format_results(results, Cpk_values, n_values, alpha_values)
    % 建立表格標題
    headers = {'Cpk', 'n', '1-α', 'Gpk', 'Bpk', 'KHpk', 'Hpk', 'Npk', ...
              'E(Gpk)', 'E(Bpk)', 'E(KHpk)', 'E(Hpk)', 'E(Npk)'};
    
    % 初始化資料矩陣
    num_rows = length(Cpk_values) * length(n_values) * length(alpha_values);
    data = zeros(num_rows, length(headers));
    row = 1;
    
    % 填充資料矩陣
    for i = 1:length(Cpk_values)
        for j = 1:length(n_values)
            for k = 1:length(alpha_values)
                temp_result = results{i}{j,k};
                data(row,:) = [Cpk_values(i), n_values(j), 1-alpha_values(k), ...
                              temp_result.coverage, temp_result.expected_values];
                row = row + 1;
            end
        end
    end
    
    % 建立並格式化表格
    formatted_table = array2table(data, 'VariableNames', headers);
    
    disp('-------------------------------------------------------------------------------------------------------------------------------');
    disp('                                     Coverage probability                           Expected value');
    disp('-------------------------------------------------------------------------------------------------------------------------------');
    disp(formatted_table);
end

simulateCpkTable()