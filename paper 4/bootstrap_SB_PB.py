import numpy as np
from scipy.stats import norm


USL = 61
LSL = 40 
T = 49
mu = 50
alpha = 0.05

sigma_list = [2, 3, 3.7]
n_list = [20, 40, 70]

B = 1000  
N = 1000  
z_value = norm.ppf(1 - alpha) 

true_cp = [(USL - LSL) / (6 * sigma) for sigma in sigma_list]
true_cpk = [min(USL - mu, mu - LSL) / (3 * sigma) for sigma in sigma_list]
true_cpm = [(USL - LSL) / (6 * (np.sqrt(sigma**2 + (mu - T)**2))) for sigma in sigma_list]

for sigma_index, sigma in enumerate(sigma_list):
    for n in n_list:
        sb_cp_count = 0
        sb_cpk_count = 0
        sb_cpm_count = 0

        pb_cp_count = 0
        pb_cpk_count = 0
        pb_cpm_count = 0

        bcpb_cp_count = 0
        bcpb_cpk_count = 0
        bcpb_cpm_count = 0



        for _ in range(N):
            sample_fixed = np.random.normal(mu, sigma, n)

            x_bar_original = np.mean(sample_fixed)
            s_original = np.std(sample_fixed, ddof=1)
            temp_original = np.sqrt(s_original**2 + (x_bar_original - T)**2)

            cp_original = (USL - LSL) / (6 * s_original)
            cpk_original = min(USL - x_bar_original, x_bar_original - LSL) / (3 * s_original)
            cpm_original = (USL - LSL) / (6 * temp_original)
            
            cp_hat_bootstrap = []
            cpk_hat_bootstrap = []
            cpm_hat_bootstrap = []

            for _ in range(B):
                sample = np.random.choice(sample_fixed, size=n, replace=True)
                x_bar = np.mean(sample)
                s = np.std(sample, ddof=1)
                temp = np.sqrt(s**2 + (x_bar - T)**2)

                cp_hat_bootstrap.append((USL - LSL) / (6 * s))
                cpk_hat_bootstrap.append(min(USL - x_bar, x_bar - LSL) / (3 * s))
                cpm_hat_bootstrap.append((USL - LSL) / (6 * temp))

            cp_hat_bootstrap.sort()
            cpk_hat_bootstrap.sort()
            cpm_hat_bootstrap.sort()

            cp_hat_bar = np.mean(cp_hat_bootstrap)
            cp_hat_s = np.std(cp_hat_bootstrap, ddof=1)
            cpk_hat_bar = np.mean(cpk_hat_bootstrap)
            cpk_hat_s = np.std(cpk_hat_bootstrap, ddof=1)
            cpm_hat_bar = np.mean(cpm_hat_bootstrap)
            cpm_hat_s = np.std(cpm_hat_bootstrap, ddof=1)

            SB_L_cp = cp_original - z_value * cp_hat_s
            SB_L_cpk = cpk_original - z_value * cpk_hat_s
            SB_L_cpm = cpm_original - z_value * cpm_hat_s

            PB_L_cp = np.percentile(cp_hat_bootstrap, alpha * 100)
            PB_L_cpk = np.percentile(cpk_hat_bootstrap, alpha * 100)
            PB_L_cpm = np.percentile(cpm_hat_bootstrap, alpha * 100)

            P0_cp = np.sum(cp_hat_bootstrap <= cp_original) / B
            P0_cpk = np.sum(cpk_hat_bootstrap <= cpk_original) / B
            P0_cpm = np.sum(cpm_hat_bootstrap <= cpm_original) / B
            
            Z0_cp = norm.ppf(P0_cp)
            Z0_cpk = norm.ppf(P0_cpk)
            Z0_cpm = norm.ppf(P0_cpm)

            PL_cp = norm.cdf(2 * Z0_cp - z_value)
            PL_cpk = norm.cdf(2 * Z0_cpk - z_value)
            PL_cpm = norm.cdf(2 * Z0_cpm - z_value)

            lower_bound_index_cp = int( PL_cp* B)
            lower_bound_index_cpk = int( PL_cpk* B)
            lower_bound_index_cpm = int( PL_cpm* B)

            BCPB_L_cp = cp_hat_bootstrap[lower_bound_index_cp]
            BCPB_L_cpk = cpk_hat_bootstrap[lower_bound_index_cpk]
            BCPB_L_cpm = cpm_hat_bootstrap[lower_bound_index_cpm]            



            sb_cp_count += SB_L_cp < true_cp[sigma_index]
            sb_cpk_count += SB_L_cpk < true_cpk[sigma_index]
            sb_cpm_count += SB_L_cpm < true_cpm[sigma_index]

            pb_cp_count += PB_L_cp < true_cp[sigma_index]
            pb_cpk_count += PB_L_cpk < true_cpk[sigma_index]
            pb_cpm_count += PB_L_cpm < true_cpm[sigma_index]

            bcpb_cp_count += BCPB_L_cp < true_cp[sigma_index]
            bcpb_cpk_count += BCPB_L_cpk < true_cpk[sigma_index]
            bcpb_cpm_count += BCPB_L_cpm < true_cpm[sigma_index]



        coverage_probability_SB_CP = sb_cp_count / N
        coverage_probability_SB_CPM = sb_cpm_count / N
        coverage_probability_SB_CPK = sb_cpk_count / N

        coverage_probability_PB_CPK = pb_cpk_count / N
        coverage_probability_PB_CP = pb_cp_count / N
        coverage_probability_PB_CPM = pb_cpm_count / N

        coverage_probability_BCPB_CPK = bcpb_cpk_count / N
        coverage_probability_BCPB_CP = bcpb_cp_count / N
        coverage_probability_BCPB_CPM = bcpb_cpm_count / N



        print(f"n: {n}, sigma: {sigma}")
        print(f"SB Cp Coverage Probability: {coverage_probability_SB_CP:.4f}")
        print(f"PB Cp Coverage Probability: {coverage_probability_PB_CP:.4f}")
        print(f"BCPB Cp Coverage Probability: {coverage_probability_BCPB_CP:.4f}")
        print()
        print(f"SB Cpk Coverage Probability: {coverage_probability_SB_CPK:.4f}")
        print(f"PB Cpk Coverage Probability: {coverage_probability_PB_CPK:.4f}")
        print(f"BCPB Cpk Coverage Probability: {coverage_probability_BCPB_CPK:.4f}")
        print()
        print(f"SB Cpm Coverage Probability: {coverage_probability_SB_CPM:.4f}")
        print(f"PB Cpm Coverage Probability: {coverage_probability_PB_CPM:.4f}")
        print(f"BCPB Cpm Coverage Probability: {coverage_probability_BCPB_CPM:.4f}")
        print()
        print()