import numpy as np
import scipy.stats as stats


# --- rMGMS fitting class ---
class Lin_rMGMS:
    def __init__(self, a_lin=1.10, b_lin=-1.95, sigma=0.20):
        self.a_lin = a_lin
        self.b_lin = b_lin
        self.sigma = sigma

    def direct(self, xi):
        return 10**self.b_lin * xi**self.a_lin
    
    def inverse(self, yi):
        return (yi / (10**self.b_lin))**(1 / self.a_lin)

# --- Statistical tests class ---    
class Statistics:
    def statistical_tests(x, y):
        # Pearson correlation
        pearson_r, pearson_p = stats.pearsonr(x, y)
        print(f"  Pearson: r = {pearson_r:.3f}, p = {pearson_p:.3f}")
        # Spearman rank correlation
        spearman_rho, spearman_p = stats.spearmanr(x, y)
        print(f"  Spearman: ρ = {spearman_rho:.3f}, p = {spearman_p:.3f}")
        # Kendall’s tau
        kendall_tau, kendall_p = stats.kendalltau(x, y)
        print(f"  Kendall: τ = {kendall_tau:.3f}, p = {kendall_p:.3f}")

        return pearson_r, pearson_p, spearman_rho, spearman_p, kendall_tau, kendall_p