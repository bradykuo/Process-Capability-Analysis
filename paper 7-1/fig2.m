% 設定基本參數
xi = 1;
cpk_values = 1:0.01:2;
m_values = [1, 2, 3, 5];
alpha = 0.05;
beta = 0.05;
CAQL = 1.33;
CLTPD = 1.00;

% 初始化ASN矩陣
ASN = zeros(length(m_values), length(cpk_values));

% 對每個m值計算n和c0，並計算ASN曲線
for i = 1:length(m_values)
    m = m_values(i);
    
    % 計算規格界限係數
    bA = 3 * CAQL + xi;
    bL = 3 * CLTPD + xi;
    
    % 定義積分函數
    S1 = @(n, c0) (1 - integral(@(t) chi2cdf((n-1)*(bA*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
        (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, bA*sqrt(n)))^m - (alpha);
    S2 = @(n, c0) (1 - integral(@(t) chi2cdf((n-1)*(bL*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
        (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, bL*sqrt(n)))^m - (1-beta);
    
    % 建立聯立方程組
    system_of_equations = @(x) [S1(x(1), x(2)); S2(x(1), x(2))];
    
    % 求解n和c0
    initial_guess = [10, 1];
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8, 'MaxIterations', 1000);
    solution = fsolve(system_of_equations, initial_guess, options);
    n = ceil(solution(1));
    c0 = solution(2);
    
    % 計算每個Cpk值的ASN
    for j = 1:length(cpk_values)
        % 計算接受機率
        b = 3*cpk_values(j) + abs(xi);
        Pa = integral(@(t) chi2cdf((n-1)*(b*sqrt(n)-t).^2/(9*n*c0^2), n-1) .* ...
            (normpdf(t + xi*sqrt(n)) + normpdf(t - xi*sqrt(n))), 0, b*sqrt(n));
        % 計算ASN
        if m == 1
            ASN(i,j) = n;
        else
            ASN(i,j) = n * (1 - (1-Pa)^m) / Pa;
        end
    end
end

% 繪圖
figure;
hold on;
styles = {'-.b', ':r', '-k', '--g'};

for i = 1:length(m_values)
    plot(cpk_values, ASN(i,:), styles{i}, 'LineWidth', 1.5);
end

% 設定圖形屬性
grid on;
xlabel('Cpk Value');
ylabel('Average Sample Number (ASN)');
axis([1 2 0 250]);
legend(['m = ' num2str(m_values(1))], ['m = ' num2str(m_values(2))], ...
       ['m = ' num2str(m_values(3))], ['m = ' num2str(m_values(4))], ...
       'Location', 'northeast');
set(gca, 'FontSize', 12);
box on;