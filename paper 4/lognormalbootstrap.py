import numpy as np
from scipy.stats import norm, lognorm, chi2

USL = 61
LSL = 40 
T = 49
mu = 50
alpha = 0.05

sigma_list = [2, 3, 3.7]
n_list = [20, 40, 70]

B = 1000  # Number of bootstrap samples
N = 1000  # Number of simulation runs
z_value = norm.ppf(1 - alpha)

# Calculate true indices
true_cp = [(USL - LSL) / (6 * sigma) for sigma in sigma_list]
true_cpk = [min(USL - mu, mu - LSL) / (3 * sigma) for sigma in sigma_list]
true_cpm = [(USL - LSL) / (6 * (np.sqrt(sigma**2 + (mu - T)**2))) for sigma in sigma_list]

def run_simulation(distribution_type='normal'):
    for sigma_index, sigma in enumerate(sigma_list):
        for n in n_list:
            sb_cp_count = sb_cpk_count = sb_cpm_count = 0
            pb_cp_count = pb_cpk_count = pb_cpm_count = 0
            bcpb_cp_count = bcpb_cpk_count = bcpb_cpm_count = 0

            for _ in range(N):
                # Generate samples based on distribution type
                if distribution_type == 'normal':
                    sample_fixed = np.random.normal(mu, sigma, n)
                elif distribution_type == 'lognormal':
                    # Convert normal parameters to lognormal parameters
                    m = np.log(mu**2 / np.sqrt(mu**2 + sigma**2))
                    s = np.sqrt(np.log(1 + (sigma**2 / mu**2)))
                    sample_fixed = np.random.lognormal(m, s, n)
                elif distribution_type == 'chisquare':
                    # Scale and shift chi-square to match desired mean and variance
                    df = 4  # degrees of freedom
                    scale = sigma / np.sqrt(2 * df)
                    shift = mu - df * scale
                    sample_fixed = chi2.rvs(df, size=n) * scale + shift

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

                # Standard Bootstrap (SB)
                cp_hat_s = np.std(cp_hat_bootstrap, ddof=1)
                cpk_hat_s = np.std(cpk_hat_bootstrap, ddof=1)
                cpm_hat_s = np.std(cpm_hat_bootstrap, ddof=1)

                SB_L_cp = cp_original - z_value * cp_hat_s
                SB_L_cpk = cpk_original - z_value * cpk_hat_s
                SB_L_cpm = cpm_original - z_value * cpm_hat_s

                # Percentile Bootstrap (PB)
                PB_L_cp = np.percentile(cp_hat_bootstrap, alpha * 100)
                PB_L_cpk = np.percentile(cpk_hat_bootstrap, alpha * 100)
                PB_L_cpm = np.percentile(cpm_hat_bootstrap, alpha * 100)

                # Bias Corrected Percentile Bootstrap (BCPB)
                P0_cp = np.sum(cp_hat_bootstrap <= cp_original) / B
                P0_cpk = np.sum(cpk_hat_bootstrap <= cpk_original) / B
                P0_cpm = np.sum(cpm_hat_bootstrap <= cpm_original) / B
                
                Z0_cp = norm.ppf(P0_cp)
                Z0_cpk = norm.ppf(P0_cpk)
                Z0_cpm = norm.ppf(P0_cpm)

                PL_cp = norm.cdf(2 * Z0_cp - z_value)
                PL_cpk = norm.cdf(2 * Z0_cpk - z_value)
                PL_cpm = norm.cdf(2 * Z0_cpm - z_value)

                lower_bound_index_cp = int(PL_cp * B)
                lower_bound_index_cpk = int(PL_cpk * B)
                lower_bound_index_cpm = int(PL_cpm * B)

                BCPB_L_cp = cp_hat_bootstrap[lower_bound_index_cp]
                BCPB_L_cpk = cpk_hat_bootstrap[lower_bound_index_cpk]
                BCPB_L_cpm = cpm_hat_bootstrap[lower_bound_index_cpm]

                # Count coverage
                sb_cp_count += SB_L_cp < true_cp[sigma_index]
                sb_cpk_count += SB_L_cpk < true_cpk[sigma_index]
                sb_cpm_count += SB_L_cpm < true_cpm[sigma_index]

                pb_cp_count += PB_L_cp < true_cp[sigma_index]
                pb_cpk_count += PB_L_cpk < true_cpk[sigma_index]
                pb_cpm_count += PB_L_cpm < true_cpm[sigma_index]

                bcpb_cp_count += BCPB_L_cp < true_cp[sigma_index]
                bcpb_cpk_count += BCPB_L_cpk < true_cpk[sigma_index]
                bcpb_cpm_count += BCPB_L_cpm < true_cpm[sigma_index]

            # Calculate coverage probabilities
            print(f"\nDistribution: {distribution_type}")
            print(f"n: {n}, sigma: {sigma}")
            print(f"SB Cp Coverage Probability: {sb_cp_count/N:.4f}")
            print(f"PB Cp Coverage Probability: {pb_cp_count/N:.4f}")
            print(f"BCPB Cp Coverage Probability: {bcpb_cp_count/N:.4f}")
            print()
            print(f"SB Cpk Coverage Probability: {sb_cpk_count/N:.4f}")
            print(f"PB Cpk Coverage Probability: {pb_cpk_count/N:.4f}")
            print(f"BCPB Cpk Coverage Probability: {bcpb_cpk_count/N:.4f}")
            print()
            print(f"SB Cpm Coverage Probability: {sb_cpm_count/N:.4f}")
            print(f"PB Cpm Coverage Probability: {pb_cpm_count/N:.4f}")
            print(f"BCPB Cpm Coverage Probability: {bcpb_cpm_count/N:.4f}")
            print("-" * 50)

# Run simulations for each distribution

print("\nRunning Lognormal distribution simulation...")
run_simulation('lognormal')
