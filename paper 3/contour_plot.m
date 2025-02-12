n_values = 10:1:150;        % 樣本數範圍
c0_values = 0.8:0.01:1.5;   % 臨界值範圍

% 初始化計算矩陣
S1_values = zeros(length(c0_values), length(n_values));
S2_values = zeros(length(c0_values), length(n_values));

alpha = 0.05;      % 生產者風險
beta = 0.1;       % 消費者風險
CAQL = 1.33;       % 可接受品質水準
CLTPD = 1.0;       % 可容忍品質水準
xi = 1;            % 製程偏移量
b1 = 3 * CAQL + xi;
b2 = 3 * CLTPD + xi;

% 計算 S1 和 S2 矩陣
for i = 1:length(c0_values)
 c0 = c0_values(i);
for j = 1:length(n_values)
 n = n_values(j);
 
 % 定義 S1 積分式
 S1 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
 .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b1 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1 - alpha);
 
 % 定義 S2 積分式
 S2 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
 .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b2 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;
 
 S1_values(i, j) = S1(n, c0);
 S2_values(i, j) = S2(n, c0);
end
end

% 建立網格
[N, C0] = meshgrid(n_values, c0_values);

% 繪製 S1 等高線圖
figure;
contour_levels_S1 = -0.9:0.1:0;
contour(N, C0, S1_values, contour_levels_S1, 'ShowText', 'on');
xlabel('n');
ylabel('C_0');
axis([10 150 0.8 1.5]);
box on;
grid off;

% 繪製 S2 等高線圖
figure;
contour_levels_S2 = 0:0.1:0.8;
contour(N, C0, S2_values, contour_levels_S2, 'ShowText', 'on');
xlabel('n');
ylabel('C_0');
axis([10 150 0.8 1.5]);
box on;
grid off;

% 求解 S1 = S2 = 0 的交點
obj_fun = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
x0 = [50, 1.2];
options = optimoptions('fsolve', 'Display', 'off');
[x_sol, fval] = fsolve(obj_fun, x0, options);
n_intersect = x_sol(1);
c0_intersect = x_sol(2);

% 繪製合併等高線圖
figure;
hold on;

% 繪製所有等高線（虛線）
[C1,h1] = contour(N, C0, S1_values, contour_levels_S1, '--k');
[C2,h2] = contour(N, C0, S2_values, contour_levels_S2, '--k');

% 繪製 S1 = 0 和 S2 = 0 的特殊線
[C3,h3] = contour(N, C0, S1_values, [0 0], 'b-', 'LineWidth', 2);
[C4,h4] = contour(N, C0, S2_values, [0 0], 'r-', 'LineWidth', 2);

% 標示交點
plot(n_intersect, c0_intersect, 'k.', 'MarkerSize', 15);
text(n_intersect+5, c0_intersect, sprintf('(n, C_0) = (%.0f, %.4f)', n_intersect, c0_intersect), 'FontSize', 10);
xlabel('n');
ylabel('C_0');
axis([10 150 0.8 1.5]);
box on;
grid off;
hold off;