% 定義範圍（反向排序）
n_values = 150:-1:10;         % n 的範圍由大到小
c0_values = 1.5:-0.01:0.8;    % c0 的範圍由大到小

% 初始化矩陣
S1_values = zeros(length(c0_values), length(n_values));
S2_values = zeros(length(c0_values), length(n_values));

% 定義參數
alpha = 0.05;
beta = 0.1;  
CAQL = 1.33;
CLTPD = 1.0;
xi = 1;

% 計算 b1 和 b2
b1 = 3 * CAQL + xi;
b2 = 3 * CLTPD + xi;

% 計算每個 (n, c0) 組合的 S1 和 S2
for i = 1:length(c0_values)
    c0 = c0_values(i);
    for j = 1:length(n_values)
        n = n_values(j);
        
        % S1 積分
        S1 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b1 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
            .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b1 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - (1 - alpha);
        
        % S2 積分
        S2 = @(n, c0) integral(@(t) chi2cdf((n - 1) * (b2 * sqrt(n) - t).^2 / (9 * n * c0^2), n - 1) ...
            .* (normpdf(t + xi * sqrt(n)) + normpdf(t - xi * sqrt(n))), 0, b2 * sqrt(n), 'RelTol', 1e-6, 'AbsTol', 1e-9) - beta;
        
        S1_values(i, j) = S1(n, c0);
        S2_values(i, j) = S2(n, c0);
    end
end

[N, C0] = meshgrid(n_values, c0_values);

% Fig 1: S1 surface plot
figure;
surf(N, C0, S1_values, 'EdgeColor', 'none');
xlabel('n');
ylabel('c0');
zlabel('S1(n, c0)');
title('Surface Plot of S1(n, c0)');
colorbar;
view(45, 30);  % 調整視角以更好地顯示曲面

% Fig 2: S2 surface plot
figure;
surf(N, C0, S2_values, 'EdgeColor', 'none');
xlabel('n');
ylabel('c0');
zlabel('S2(n, c0)');
title('Surface Plot of S2(n, c0)');
colorbar;
view(45, 30);  % 調整視角以更好地顯示曲面

% Fig 3: Combined surface plot
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
view(45, 30);  % 調整視角以更好地顯示曲面
hold off;