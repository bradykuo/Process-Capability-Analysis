% 計算給定 Cpk 值的 PPM 界限
function [lower_ppm, upper_ppm] = calculate_ppm_bounds(cpk)
% 單邊規格計算
 upper_ppm = (1 - normcdf(3 * cpk)) * 1e6;
% 雙邊規格計算（製程中心對稱情況）
 lower_ppm = 2 * (1 - normcdf(3 * cpk)) * 1e6;
end

cpk_values = [0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.24, 1.25, ...
 1.30, 1.33, 1.40, 1.45, 1.50, 1.60, 1.67, 1.70, 1.80, 1.90, 2.00];

% 計算每個 Cpk 值的 PPM 界限
results = zeros(length(cpk_values), 2);
for i = 1:length(cpk_values)
 [results(i,1), results(i,2)] = calculate_ppm_bounds(cpk_values(i));
end

results = round(results, 3);

% 建立並顯示表格
T = table(cpk_values', results(:,1), results(:,2), ...
'VariableNames', {'Cpk', 'Upper_Bound', 'Lower_Bound'});
disp(T)

% 驗證特定值與論文的對照
fprintf('\nVerification of selected values:\n');
[lower_1_33, upper_1_33] = calculate_ppm_bounds(1.33);
fprintf('For Cpk = 1.33: Expected ≈ 66/33 PPM, Calculated: %.3f %.3f PPM\n', ...
 lower_1_33, upper_1_33);
[lower_1_67, upper_1_67] = calculate_ppm_bounds(1.67);
fprintf('For Cpk = 1.67: Expected ≈ 0.544/0.272 PPM, Calculated: %.3f %.3f PPM\n', ...
 lower_1_67, upper_1_67);