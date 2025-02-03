% 設定參數
xi = 1;
cpk_values = 1:0.01:2;
m_values = [1, 2, 3, 5];
n = 100;    % 固定 n 值
c0 = 1.5;  % 固定 c0 值

% 初始化
PA = zeros(length(m_values), length(cpk_values));

% 接受機率計算函數
function Pa = calc_acceptance_prob(n, c0, cpk, xi)
    b = 3*cpk + abs(xi);
    Pa = integral(@(t) chi2cdf((n-1)*(b*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
        (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, b*sqrt(n));
end

% 計算每個m值的接受機率
for i = 1:length(m_values)
    m = m_values(i);
    
    % 計算接受機率
    for j = 1:length(cpk_values)
        Pa = calc_acceptance_prob(n, c0, cpk_values(j), xi);
        PA(i,j) = 1 - (1-Pa)^m;
    end
end

% 繪圖
figure;
hold on;
styles = {'-.b', ':r', '-k', '--g'};
for i = 1:length(m_values)
    plot(cpk_values, PA(i,:), styles{i}, 'LineWidth', 1.5);
end
grid on;
xlabel('Cpk Value');
ylabel('Probability of acceptance');
axis([1 2 0 1]);
legend(['m = ' num2str(m_values(1))], ['m = ' num2str(m_values(2))], ...
    ['m = ' num2str(m_values(3))], ['m = ' num2str(m_values(4))], ...
    'Location', 'southeast');
set(gca, 'FontSize', 12);
box on;