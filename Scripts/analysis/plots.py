import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np
import matplotlib
from matplotlib.ticker import LogLocator, LogFormatterExponent
from cycler import cycler


# --- MNRAS plot formatting ---
def set_rc_params(mult):
    matplotlib.rcParams.update({'font.size': 20*mult})
    matplotlib.rcParams['legend.fontsize'] = 17.5*mult
    matplotlib.rcParams['axes.linewidth'] = 1.5
    matplotlib.rcParams['xtick.labelsize'] = 20*mult
    matplotlib.rcParams['ytick.labelsize'] = 20*mult
    matplotlib.rcParams['xtick.major.size'] = 10
    matplotlib.rcParams['ytick.major.size'] = 10
    matplotlib.rcParams['xtick.major.width'] = 2
    matplotlib.rcParams['ytick.major.width'] = 2
    matplotlib.rcParams['xtick.minor.size'] = 5
    matplotlib.rcParams['ytick.minor.size'] = 5
    matplotlib.rcParams['xtick.minor.width'] = 1
    matplotlib.rcParams['ytick.minor.width'] = 1
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.bottom'] = True
    matplotlib.rcParams['ytick.left'] = True
    matplotlib.rcParams['axes.labelsize'] = 20*mult

# --- Class containing plotting methods ---
class Plotter:
    def __init__(self, show_figures):
        self.show_figures = show_figures
        set_rc_params(0.5) # setting MNRAS style with a multiplier of 0.8 for font sizes
        
    def fits_image(self, data, title):
        plt.imshow(data, origin='lower', cmap='gray')
        plt.colorbar()
        plt.title(title)
        if self.show_figures:
            plt.show()
        plt.close()

    # Make plot using two lists from a dictionary of variables
    def dic_plot(self, dic_plot, variable_1, variable_2, variable_2_err):
        fig, ax = plt.subplots(figsize=(6,6))
        sc = ax.scatter(dic_plot[variable_1][1:], dic_plot[variable_2][1:], s=5.0, alpha=1.0, color = 'black') #, label="Gradient vs Resolution") # Exclude element 0 (units)
        # Scatter with error bars
        ax.errorbar(
            dic_plot[variable_1][1:], dic_plot[variable_2][1:], yerr=dic_plot[variable_2_err][1:],
            fmt='o', ms=2.5, color='black', ecolor='gray', elinewidth=0.8, capsize=2.0, alpha=0.8
        )
        slope, intercept, r_value, p_value, std_err = linregress(dic_plot[variable_1][1:], dic_plot[variable_2][1:])
        x_fit = np.linspace(np.min(dic_plot[variable_1][1:]), np.max(dic_plot[variable_1][1:]), 200)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, color='red', lw=1.5, label=f'Fit: y={slope:.5f}x+{intercept:.2f}\n$R^2$={r_value**2:.2f}')
        ax.set_xlabel(f'{variable_1} [{dic_plot[variable_1][0]}]')
        ax.set_ylabel(f'{variable_2} [{dic_plot[variable_2][0]}]')
        ax.legend()
        fig.tight_layout()
        figure = f'./CompareRuns/{variable_1}_{variable_2}'
        fig.savefig(figure, dpi=300)
        if self.show_figures:
            plt.show()
        plt.close()

    def plot_parameters(self, variable_1, variable_1_name, variable_1_units, variable_2, variable_2_name, variable_2_units, path):    
        fig, ax = plt.subplots()
        # Scatter plot
        ax.scatter(variable_1, variable_2, c='black', marker='.', s=5.0, alpha=1.0)
        # Best-fit linear regression
        slope, intercept, r_value, p_value, std_err = linregress(variable_1, variable_2)
        x_fit = np.linspace(np.min(variable_1), np.max(variable_1), 200)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, color='red', lw=1.5, label=f'Fit: y={slope:.2f}x+{intercept:.2f}\n$R^2$={r_value**2:.2f}')
        ax.set_xlabel(f'{variable_1_name} [{variable_1_units}]') # omit log since it is already in dex and it's an offset (it's implied)
        ax.set_ylabel(f'{variable_2_name} [{variable_2_units}]')
        ax.legend(loc='best', fontsize='small')
        fig.tight_layout()
        rMGMS_figure = f'{path}/WPX_Figures/{variable_1_name}_vs_{variable_2_name}.png'
        fig.savefig(rMGMS_figure, dpi=300)
        if self.show_figures:
            plt.show()
        plt.close()

    def plot_rMGMS(self, a, b, sig, ro, path, x, y, x_pivot):
        # Plane defined by Lin et al. (2020): log(Σ_H2) = alpha*log(Σ_*) + beta*log(Σ_SFR) + C
        alpha = 1.10
        beta = 1/1.05
        C = 15.23
        # The projection onto the rMGMS (Σ_H2 vs. Σ_*) is:
        # log(Σ_H2) = a*log(Σ_*) + beta*log(Σ_SFR) + C
        # Below are various fixed values for log(Σ_SFR) in order to plot several lines of the rMGMS projection from the plane
        log_SFR_values = [-21, -20, -19, -18, -17, -16]
        # Generate values according to x and y
        x_lin = np.logspace(9, 12.1, 200)
        y_lin = 10**b * x_lin**a
        sig_mult = 1
        y_lin_upper = y_lin * 10 ** (sig_mult*sig)
        y_lin_lower = y_lin * 10 ** (sig_mult*-sig)
        fig, ax = plt.subplots()
        # set log scales
        ax.set_xscale('log')
        ax.set_yscale('log')
        # major ticks at powers of ten
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_locator(LogLocator(base=10))
        # formatter that only shows the exponent
        ax.xaxis.set_major_formatter(LogFormatterExponent(base=10))
        ax.yaxis.set_major_formatter(LogFormatterExponent(base=10))
        # Plot Lin's fit
        ax.plot(x_lin, y_lin, color = 'black', lw = 1, label=f'Lin et al. (2020): a = {a}; b = {b}')
        ax.plot(x_lin, y_lin_upper, color = 'black', linestyle = '--', lw = 0.7, label = fr'$+\sigma$')
        ax.plot(x_lin, y_lin_lower, color = 'black', linestyle = '--', lw = 0.7, label = fr'$-\sigma$')
        # Check if x and y are lists (all galaxies together) or dictionaries (each galaxy separately)
        if isinstance(y, np.ndarray):
            ax.scatter(x, y, c='black', marker='.', s=10, alpha=1.0, linewidths=0)
            # --- Best-fit linear regression ---
            # Fit a line in log–log space
            logx = np.log10(x)  
            logy = np.log10(y)
            # Pivoting = 0 or mean of x.
            # Small changes in m don’t require large changes in c to preserve the same fit near your data.
            # Parameter correlation drops sharply.
            # The regression is numerically better conditioned (smaller rounding errors, more meaningful error bars).
            logx_minus_x_pivot = logx if x_pivot == 0 else logx - np.log10(x_pivot)
            slope, intercept, r_value, p_value, std_err_slope = linregress(logx_minus_x_pivot, logy)
            n = len(logx_minus_x_pivot)
            x_mean = np.mean(logx_minus_x_pivot)
            s_yx = std_err_slope * np.sqrt(np.sum((logx_minus_x_pivot - x_mean)**2))  # equivalent to s
            std_err_intercept = s_yx * np.sqrt(1/n + x_mean**2 / np.sum((logx_minus_x_pivot - x_mean)**2))
            # Extend line from 10^5 to 10^12 in x
            x_fit = np.logspace(9, 12.1, 200)
            # power-law fit in linear space, divide by x_pivot to include the pivot correction
            y_fit = 10**intercept * x_fit**slope if x_pivot == 0 else 10**intercept * (x_fit / x_pivot)**slope
            # --- Actual plotting ---
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.plot(x_fit, y_fit, color='red', lw=1.5,
                label=f'Fit: log(y)={slope:.2f}[log(x)-log(x_pivot)]+{intercept:.2f}\n$R^2$={r_value**2:.2f}')
            for log_SFR in log_SFR_values:
                y_plane = 10**(C + beta*log_SFR) * x_fit**alpha
                ax.plot(x_fit, y_plane, linestyle=':', lw=1,
                        label=f'Plane projection (log Σ_SFR={log_SFR})')
                        #label=f'Plane projection (log Σ_SFR={log_SFR}): a={alpha}, b={beta:.2f}')
            ax.set_xlabel(r'Log $\Sigma_{*}$ [M$_\odot$ kpc$^{-2}$]')
            ax.set_ylabel(r'Log $\Sigma_{\rm H_2}$ [M$_\odot$ kpc$^{-2}$]')
            #ax.set_title('Pixel-by-Pixel: Σ_H2 vs. Σ_*')
            ax.legend()
            fig.tight_layout()
            rMGMS_figure = f'{path}/rMGMS_Figures/Total_rMGMS_vs_Lin2020_log_monochrome_pivot{x_pivot:.3e}.png'
            fig.savefig(rMGMS_figure, dpi=300)
            if self.show_figures:
                plt.show()
            plt.close()
            
            return slope, intercept, std_err_slope, std_err_intercept

        # Each galaxy is a different color
        else:
            ax.set_prop_cycle(cycler(color=plt.cm.Set1.colors))
            for key in x:
                sig_star = np.genfromtxt(x[key][1:], usecols=0, dtype=float) # Sigma_star
                sig_h2 = np.genfromtxt(x[key][1:], usecols=1, dtype=float) # Sigma_H2
                # Linewidths for dot-like effect, '.' is not enough. Good for clouds. Only for the png file
                ax.scatter(sig_star, sig_h2, marker='.', s=10, alpha=1.0, linewidths=0, label=f'{key}')
                #ax.scatter(sig_star, sig_h2, marker='.', linewidths=0, label=f'{key}')
            ax.set_xlabel(r'Log $\Sigma_{*}$ [M$_\odot$ kpc$^{-2}$]')
            ax.set_ylabel(r'Log $\Sigma_{\rm H_2}$ [M$_\odot$ kpc$^{-2}$]')
            #ax.set_title('Pixel-by-Pixel: Σ_H2 vs. Σ_*')
            ax.legend()
            fig.tight_layout()
            rMGMS_figure = f'{path}/rMGMS_Figures/Total_rMGMS_vs_Lin2020_logm_multichrome.png'
            fig.savefig(rMGMS_figure, dpi=300)
            if self.show_figures:
                plt.show()
            plt.close()

    def plot_center(self, data, center_x, center_y):
        plt.imshow(data, origin='lower', cmap='gray')
        plt.colorbar()
        plt.scatter(center_x, center_y, color='red', marker='.', s=25, label='Galaxy Center')
        plt.title("Imported Galaxy center H2")
        if self.show_figures:
            plt.show()
        plt.close()