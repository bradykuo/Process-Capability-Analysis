function plotCpkAnalysis()
    % 開始計時
    tic;
    
    % 設定參數
    n_values = [30, 50, 70, 100, 150, 200];      % 樣本數範圍
    xi_values = 0:0.1:3;                         % xi 值範圍
    Cpk_hat_values = [0.7, 0.9, 1.2, 2.0, 2.5, 3.0];  % Cpk 估計值
    
    % 創建圖形視窗並設定大小
    figure('Position', [100, 100, 1200, 800]);
    % 創建 2x3 的子圖布局
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 為不同的 n 值設定顏色方案
    colors = lines(length(n_values));
    
    % 創建子圖
    for i = 1:length(Cpk_hat_values)
        % 計算當前 Cpk 值的所有結果
        results = calculate_C_vs_xi(Cpk_hat_values(i), n_values, xi_values);
        
        % 創建新的子圖
        nexttile
        hold on
        % 為每個樣本量繪製一條線
        for j = 1:length(n_values)
            idx = results.n == n_values(j);
            plot(results.xi(idx), results.C(idx), 'Color', colors(j,:), 'LineWidth', 1.5)
        end
        hold off
        
        % 設定子圖標題和軸標籤
        title(['Ĉpk = ', num2str(Cpk_hat_values(i))])
        xlabel('|ξ|')
        ylabel('C')
        grid on
        box on
    end
    
    % 添加標題
    title(t, 'C vs |ξ| for different Ĉpk values', 'FontSize', 16, 'FontWeight', 'bold')
    
    % 添加共同圖例
    leg = legend(string(n_values), 'Location', 'southoutside', 'Orientation', 'horizontal');
    leg.Layout.Tile = 'south';
    title(leg, 'n')
    
    % 顯示執行時間
    execution_time = toc;
    fprintf('Total execution time: %.2f seconds\n', execution_time);
end

% 計算不同 xi 值和樣本數組合的 C 值
function results = calculate_C_vs_xi(Cpk_hat, n_values, xi_values)
    % 初始化結果
    [N, Xi] = meshgrid(n_values, xi_values);
    results.n = N(:);
    results.xi = Xi(:);
    results.C = zeros(size(results.n));
    
    % 計算每種組合的 C 值
    for i = 1:length(results.n)
        results.C(i) = calculate_lcb(results.n(i), Cpk_hat, results.xi(i));
    end
end

% 計算置信下界
function lcb = calculate_lcb(n, Cpk_hat, xi_hat, gamma)
    if nargin < 4
        gamma = 0.95;    % 默認置信水平
    end
    
    % 使用 fzero 尋找根
    options = optimset('Display', 'off');
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
plotCpkAnalysis()