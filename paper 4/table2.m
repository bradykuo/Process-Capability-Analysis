% 定義常數
USL = 61;          % 規格上限
LSL = 40;          % 規格下限
T = 49;            % 目標值
alpha = 0.05;      % 顯著水準

% 定義測試參數陣列
mu_list = [50, 52];           
sigma_list = [2.0, 3.0, 3.7]; 
n_list = [20, 40, 70];        

B = 1000;  % Bootstrap 重抽樣次數
N = 1000;  % 模擬重複次數
z_value = norminv(1 - alpha); % 計算 z 值

% 主要模擬迴圈
for mu_idx = 1:length(mu_list)
    mu = mu_list(mu_idx);
    
    for sigma_idx = 1:length(sigma_list)
        sigma = sigma_list(sigma_idx);
        
        % 計算真實製程能力指標值
        true_cp = (USL - LSL)/(6 * sigma);     
        true_cpk = min((USL - mu)/(3 * sigma), (mu - LSL)/(3 * sigma)); 
        true_cpm = (USL - LSL)/(6 * sqrt(sigma^2 + (mu - T)^2));        
        
        for n_idx = 1:length(n_list)
            n = n_list(n_idx);
            
            % 初始化計數器
            sb_counts = zeros(1,3);   % Standard Bootstrap
            pb_counts = zeros(1,3);    % Percentile Bootstrap
            bcpb_counts = zeros(1,3);  % Biased Corrected Percentile Bootstrap
            
            for sim = 1:N
                % 生成原始樣本
                sample_fixed = normrnd(mu, sigma, [1,n]);
                
                % 計算原始統計量
                x_bar_original = mean(sample_fixed);
                s_original = std(sample_fixed);
                temp_original = sqrt(s_original^2 + (x_bar_original - T)^2);
                
                % 計算原始製程能力指標
                cp_original = (USL - LSL)/(6 * s_original);
                cpk_original = min((USL - x_bar_original)/(3 * s_original), ...
                                 (x_bar_original - LSL)/(3 * s_original));
                cpm_original = (USL - LSL)/(6 * temp_original);
                
                % 初始化 Bootstrap 陣列
                cp_boot = zeros(1,B);
                cpk_boot = zeros(1,B);
                cpm_boot = zeros(1,B);
                
                % Bootstrap 重抽樣迴圈
                for b = 1:B
                    boot_indices = randi(n, 1, n);  % 隨機抽樣索引
                    sample = sample_fixed(boot_indices);
                    
                    % 計算 Bootstrap 樣本統計量
                    x_bar = mean(sample);
                    s = std(sample);
                    temp = sqrt(s^2 + (x_bar - T)^2);
                    
                    % 計算 Bootstrap 製程能力指標
                    cp_boot(b) = (USL - LSL)/(6 * s);
                    cpk_boot(b) = min((USL - x_bar)/(3 * s), (x_bar - LSL)/(3 * s));
                    cpm_boot(b) = (USL - LSL)/(6 * temp);
                end
                
                % 排序 Bootstrap 估計值
                cp_boot = sort(cp_boot);
                cpk_boot = sort(cpk_boot);
                cpm_boot = sort(cpm_boot);
                
                % 計算 Standard Bootstrap 信賴區間下限
                cp_boot_std = std(cp_boot);
                cpk_boot_std = std(cpk_boot);
                cpm_boot_std = std(cpm_boot);
                
                SB_L_cp = cp_original - z_value * cp_boot_std;
                SB_L_cpk = cpk_original - z_value * cpk_boot_std;
                SB_L_cpm = cpm_original - z_value * cpm_boot_std;
                
                % 計算 Percentile Bootstrap 信賴區間下限
                PB_L_cp = prctile(cp_boot, alpha * 100);
                PB_L_cpk = prctile(cpk_boot, alpha * 100);
                PB_L_cpm = prctile(cpm_boot, alpha * 100);
                
                % 計算 BCPB
                % 計算原始估計值小於等於 Bootstrap 估計值的比例
                P0_cp = sum(cp_boot <= cp_original)/B;
                P0_cpk = sum(cpk_boot <= cpk_original)/B;
                P0_cpm = sum(cpm_boot <= cpm_original)/B;
                
                % 計算標準常態分位數
                Z0_cp = norminv(P0_cp);
                Z0_cpk = norminv(P0_cpk);
                Z0_cpm = norminv(P0_cpm);
                
                % 計算偏差校正機率
                PL_cp = normcdf(2 * Z0_cp - z_value);
                PL_cpk = normcdf(2 * Z0_cpk - z_value);
                PL_cpm = normcdf(2 * Z0_cpm - z_value);
                
                % 計算對應的 Bootstrap 樣本索引
                idx_cp = max(1, round(PL_cp * B));
                idx_cpk = max(1, round(PL_cpk * B));
                idx_cpm = max(1, round(PL_cpm * B));
                
                % 取得 BCPB 信賴區間下限
                BCPB_L_cp = cp_boot(idx_cp);
                BCPB_L_cpk = cpk_boot(idx_cpk);
                BCPB_L_cpm = cpm_boot(idx_cpm);
                
                % 更新涵蓋率計數
                sb_counts = sb_counts + [SB_L_cp < true_cp, ...
                                       SB_L_cpk < true_cpk, ...
                                       SB_L_cpm < true_cpm];
                
                pb_counts = pb_counts + [PB_L_cp < true_cp, ...
                                       PB_L_cpk < true_cpk, ...
                                       PB_L_cpm < true_cpm];
                
                bcpb_counts = bcpb_counts + [BCPB_L_cp < true_cp, ...
                                           BCPB_L_cpk < true_cpk, ...
                                           BCPB_L_cpm < true_cpm];
            end
            
            % 計算涵蓋機率
            cp_prob_SB = sb_counts/N;
            cp_prob_PB = pb_counts/N;
            cp_prob_BCPB = bcpb_counts/N;
            
            % 顯示結果
            fprintf('mu: %.1f, sigma: %.1f, n: %d\n', mu, sigma, n);
            fprintf('True Values - Cp: %.2f, Cpk: %.2f, Cpm: %.2f\n', ...
                    true_cp, true_cpk, true_cpm);
            fprintf('Coverage Probabilities:\n');
            fprintf('        Cp      Cpk     Cpm\n');
            fprintf('SB:     %.3f   %.3f   %.3f\n', cp_prob_SB);
            fprintf('PB:     %.3f   %.3f   %.3f\n', cp_prob_PB);
            fprintf('BCPB:   %.3f   %.3f   %.3f\n\n', cp_prob_BCPB);
        end
    end
end