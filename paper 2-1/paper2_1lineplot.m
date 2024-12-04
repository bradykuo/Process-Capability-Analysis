% 通過數值積分計算 c0 的函數
function c0 = calculate_c0(C, Cp, n, alpha)
% 設定根查找的初始邊界
lower = 0.5;
upper = 5;
max_attempts = 10;
attempts = 0;

% 依需要調整邊界值
while attempts < max_attempts
obj_lower = objective_func(lower, C, Cp, n, alpha);
obj_upper = objective_func(upper, C, Cp, n, alpha);
   
   % 檢查是否找到合適的區間
   if ~isnan(obj_lower) && ~isnan(obj_upper) && obj_lower * obj_upper < 0
       break;
   end
   
   % 根據函數值調整邊界
   if isnan(obj_lower) || (~isnan(obj_lower) && obj_lower > 0)
       lower = lower / 2;
   end
   if isnan(obj_upper) || (~isnan(obj_upper) && obj_upper < 0)
       upper = upper * 2;
   end
   
   attempts = attempts + 1;
end

% 使用 fzero 查找根值
try
   options = optimset('Display', 'off');
   c0 = fzero(@(c0) objective_func(c0, C, Cp, n, alpha), [lower, upper], options);
catch
   c0 = NaN;
end
end

% 數值積分的目標函數
function result = objective_func(c0, C, Cp, n, alpha)
try
   % 定義積分上限
   upper_limit = 3 * Cp * sqrt(n);
   
   % 執行數值積分
   result = integral(@(y) integrand(y, c0, C, Cp, n), 0, upper_limit) - alpha;
catch
   result = NaN;
end
end

% 被積函數
function val = integrand(y, c0, C, Cp, n)
% 計算 G 項（卡方分配的累積分配函數）
chi_arg = (n - 1) * (3 * Cp * sqrt(n) - y).^2 ./ (9 * n * c0^2);
G_term = chi2cdf(chi_arg, n - 1);

% 計算 f 項（常態分配密度函數的和）
f_term = normpdf(y + 3 * (Cp - C) * sqrt(n)) + normpdf(y - 3 * (Cp - C) * sqrt(n));

% 組合兩項
val = G_term .* f_term;
end

function plot_cpk_analysis()
% 設定 Cpk 值和樣本大小
Cpk_values = [1.00, 1.33, 1.50, 1.67, 2.00];
n_values = [30, 50, 70, 100, 150, 200, 250, 300];
alpha = 0.05;

% 創建圖形視窗
figure('Position', [100, 100, 1200, 800]);

% 為每個 Cpk 值創建圖形
for i = 1:length(Cpk_values)
   % 創建子圖
   subplot(3, 2, i);
   
   % 生成 Cp 值範圍
   Cp_values = linspace(Cpk_values(i), Cpk_values(i) + 1, 20);
   
   % 初始化 c0 值矩陣
   c0_values = zeros(length(n_values), length(Cp_values));
   
   % 計算每種組合的 c0 值
   for j = 1:length(n_values)
       for k = 1:length(Cp_values)
           c0_values(j,k) = calculate_c0(Cpk_values(i), Cp_values(k), n_values(j), alpha);
       end
   end
   
   % 繪製曲線
   hold on;
   colors = lines(length(n_values));
   for j = 1:length(n_values)
       plot(Cp_values, c0_values(j,:), 'LineWidth', 1, 'Color', colors(j,:));
   end
   hold off;
   
   % 設定圖形屬性
   title(sprintf('Cpk = %.2f', Cpk_values(i)));
   xlabel('Cp');
   ylabel('c0');
   grid on;
   
   % 調整 y 軸範圍
   ylim_curr = ylim;
   ylim([ylim_curr(1), ylim_curr(2)]);
   
   % 只在第一個圖加入圖例
   if i == 1
       legend(cellstr(num2str(n_values', 'n=%d')), 'Location', 'eastoutside');
   end
end

% 加入標題
sgtitle('Plots of c0 versus Cp for different Cpk values (α = 0.05)', 'FontSize', 14);

% 調整子圖間距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.8, 0.8]);
end

plot_cpk_analysis()  % 執行繪圖函數