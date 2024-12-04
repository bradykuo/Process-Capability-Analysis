% 定義製程參數
USL = 61; % 規格上限
LSL = 40; % 規格下限
T = 49;   % 目標值

% 定義不同的平均值和標準差組合
mu = [50 52 50 52 50 52];     
sigma = [2.0 2.0 3.0 3.0 3.7 3.7];

% 初始化製程能力指標陣列
Cp = zeros(size(mu));  
Cpk = zeros(size(mu));  
Cpm = zeros(size(mu));  

% 計算每種組合的製程能力指標
for i = 1:length(mu)
    % 計算Cp (製程能力指數)
    Cp(i) = (USL - LSL)/(6*sigma(i));
    
    % 計算Cpk (考慮製程平均值偏移的製程能力指數)
    Cpu = (USL - mu(i))/(3*sigma(i));  % 上限製程能力
    Cpl = (mu(i) - LSL)/(3*sigma(i));  % 下限製程能力
    Cpk(i) = min(Cpu, Cpl);
    
    % 計算Cpm (考慮目標值的 Taguchi 製程能力指數)
    sigma_squared = sigma(i)^2 + (mu(i) - T)^2;  % 變異數加上偏離目標的平方
    Cpm(i) = (USL - LSL)/(6*sqrt(sigma_squared));
end

% 建立結果表格
data = [mu' sigma' Cp' Cpk' Cpm'];
columnNames = {'μ', 'σ', 'Cp', 'Cpk', 'Cpm'};
T = array2table(data, 'VariableNames', columnNames);

% 顯示表格
disp(T)