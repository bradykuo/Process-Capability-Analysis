function plot_cpk_analysis()
    % 創建主視窗
    figure('Position', [100, 100, 1200, 800]);

    % 設定參數
    C_values = [1.00, 1.33, 1.50, 1.67, 2.00];  % Cpk值
    Cp_ranges = {[1.00, 2.00],
                 [1.33, 2.33],
                 [1.50, 2.50],
                 [1.67, 2.67],
                 [2.00, 3.00]};
    alpha = 0.05;

    % 為每個 Cpk 值創建子圖
    for i = 1:length(C_values)
        % 創建網格數據
        n_values = linspace(30, 300, 50);
        Cp_values = linspace(Cp_ranges{i}(1), Cp_ranges{i}(2), 50);
        [N, CP] = meshgrid(n_values, Cp_values);
        
        % 計算 c0 值
        c0_values = zeros(size(N));
        for j = 1:size(N, 1)
            for k = 1:size(N, 2)
                c0_values(j,k) = calculate_c0(C_values(i), CP(j,k), N(j,k), alpha);
            end
        end

        % 創建子圖
        subplot(2, 3, i);
        surf(N, CP, c0_values);
        
        % 設定圖形屬性
        title(sprintf('Cpk = %.2f', C_values(i)));
        xlabel('n');
        ylabel('Cp');
        zlabel('c0');
        colorbar;
        view(45, 30);
        grid on;
    end

    % 添加總標題
    sgtitle('Surface plots of c0 for different Cpk values (α = 0.05)', 'FontSize', 14);
end

function c0 = calculate_c0(C, Cp, n, alpha)
    % 設定初始搜索區間
    lower = 0.5;
    upper = 5;
    
    % 調整搜索區間
    while objective_func(lower, C, Cp, n, alpha) * objective_func(upper, C, Cp, n, alpha) > 0
        if objective_func(lower, C, Cp, n, alpha) > 0
            lower = lower / 2;
        else
            upper = upper * 2;
        end
    end
    
    % 使用 fzero 找出根值
    options = optimset('Display', 'off');
    c0 = fzero(@(x) objective_func(x, C, Cp, n, alpha), [lower, upper], options);
end

function result = objective_func(c0, C, Cp, n, alpha)
    % 定義積分上限
    upper_limit = 3 * Cp * sqrt(n);
    
    % 執行數值積分
    result = integral(@(y) integrand(y, c0, C, Cp, n), 0, upper_limit) - alpha;
end

function val = integrand(y, c0, C, Cp, n)
    % 計算卡方分配項
    chi_arg = (n - 1) * (3 * Cp * sqrt(n) - y).^2 ./ (9 * n * c0^2);
    G_term = chi2cdf(chi_arg, n - 1);
    
    % 計算常態分配項
    f_term = normpdf(y + 3 * (Cp - C) * sqrt(n)) + normpdf(y - 3 * (Cp - C) * sqrt(n));
    
    % 組合結果
    val = G_term .* f_term;
end

plot_cpk_analysis()
