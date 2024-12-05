n_values = 150:-1:10;      % 樣本數範圍
c0_values = 1.5:-0.01:0.8; % 臨界值範圍

S1_values = zeros(length(c0_values), length(n_values));
S2_values = zeros(length(c0_values), length(n_values));

alpha = 0.05;     % 生產者風險
beta = 0.1;       % 消費者風險
CAQL = 1.33;      % 可接受品質水準
CLTPD = 1.0;      % 容忍品質水準
xi = 1;           % 製程偏移量
b1 = 3 * CAQL + xi;
b2 = 3 * CLTPD + xi;

% 計算 S1 和 S2 數值
for i = 1:length(c0_values)
 c0 = c0_values(i);
for j = 1:length(n_values)
 n = n_values(j);
 
 % S1 和 S2 積分計算
 S1 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
 .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b1 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1 - alpha);
 
 S2 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
 .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b2 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;
 
 S1_values(i, j) = S1(n, c0);
 S2_values(i, j) = S2(n, c0);
end
end

[N, C0] = meshgrid(n_values, c0_values);

% 繪製 S1 曲面圖
figure;
surf(N, C0, S1_values, 'EdgeColor', 'none');
xlabel('n');
ylabel('c0');
zlabel('S1(n, c0)');
title('Surface Plot of S1(n, c0)');
colorbar;
view(45, 30);

% 繪製 S2 曲面圖
figure;
surf(N, C0, S2_values, 'EdgeColor', 'none');
xlabel('n');
ylabel('c0');
zlabel('S2(n, c0)');
title('Surface Plot of S2(n, c0)');
colorbar;
view(45, 30);

% 繪製合併曲面圖
figure;
hold on;
surf(N, C0, S1_values, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'b');
surf(N, C0, S2_values, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'r');
xlabel('n');
ylabel('c0');
zlabel('S1 and S2 values');
title('3D Surface Plot of S1(n, c0) and S2(n, c0)');
legend({'S1(n, c0)', 'S2(n, c0)'}, 'Location', 'best');
colorbar;
view(45, 30);
hold off;