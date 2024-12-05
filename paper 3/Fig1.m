% 計算給定 Cpk 值的 PPM 界限
function [lower_ppm, upper_ppm] = calculate_ppm_bounds(cpk)
% 單邊規格計算
 upper_ppm = (1 - normcdf(3 * cpk)) * 1e6;
% 雙邊規格計算（製程中心對稱情況）
 lower_ppm = 2 * (1 - normcdf(3 * cpk)) * 1e6;
end

cpk_seq = 0.6:0.01:2.0;

% 計算每個 Cpk 值的 PPM 界限
results = zeros(length(cpk_seq), 2);
for i = 1:length(cpk_seq)
 [results(i,1), results(i,2)] = calculate_ppm_bounds(cpk_seq(i));
end

% 繪製圖形
figure;
plot(cpk_seq, results(:,1), 'LineWidth', 1);
hold on;
plot(cpk_seq, results(:,2), 'LineWidth', 1);

% 設定座標軸標籤和範圍
xlabel('C_{pk}', 'FontSize', 12);
ylabel('PPM', 'FontSize', 12);
ylim([0 80000]);
xlim([0.6 2.0]);

% 設定 y 軸刻度
yticks(0:20000:80000);

% 調整圖形邊界
set(gca, 'Position', [0.15 0.2 0.8 0.7]);
hold off;