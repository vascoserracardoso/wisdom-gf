#!/usr/bin/env python3

# --- Import packages ---
import os
import numpy as np
from ltsfit import ltsfit  # exact API name may vary; see package docs
import sys

# --- Import custom modules ---
from .analysis import Statistics, Plotter, Lin_rMGMS
from .Dictionary import galaxies
from .utils import userChoice

# --- Command line interface ---
def cli():
    imported_add = False
    runs_folder = './Runs'
    options = [name for name in sorted(os.listdir(runs_folder)) if os.path.isdir(os.path.join(runs_folder, name))]
    try:
        run_name = userChoice(options).runExclusiveChoice()
    except:
        print("Error -> There are no SINGLE RUNS in ./Runs. This script works only on SINGLE RUNS. It does not work on RUNSETS.")
        print("         Perform a SINGLE RUN first !")
        sys.exit()
    show_figures = input("\nShow figures? (y/n): ").strip().lower()
    if show_figures == "y":
        show_figures = True
    else:
        show_figures = False
    main(imported_add,
         run_name,
         show_figures)


def main(imported_script, run_name, show_figures):

    # Initialize Plotter
    plotter = Plotter(show_figures)

    # Parameters for the best fit in Lin et al. (2020)
    a = 1.10
    b = -1.95
    sig = 0.20 
    ro = 0.76 
    print(f'\n- Using Lin et al. (2020) fit: a = {a:.2f}, b = {b:.2f}, σ = {sig:.2f}, ρ = {ro:.2f}')

    # Initialize Lin_rMGMS class
    lin = Lin_rMGMS()

    # --- Handle directories and paths ---
    newpath = f'./Runs/{run_name}/runOutputs' 
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        os.makedirs(f'{newpath}/WPX_Figures')
        os.makedirs(f'{newpath}/rMGMS_Figures')
    base_folder = './Runs/' + run_name + '/Galaxies'
    filename = 'rMGMSOutputs.txt'
    output_file = f'{newpath}/All_rMGMSOutputs.txt'

    # --- Compile all rMGMSOutputs.txt into a single file and means of Sigma_star and Sigma_H2 for each galaxy ---
    all_lines = []
    header_written = False
    rMGMS_galaxy_dic = {}
    print("- Offsets in relation to Lin+2020 rMGMS...")
    galaxies_offset_star_dex = {}
    galaxies_offset_h2_dex = {}
    for dirname in sorted(os.listdir(base_folder)):
        # One galaxy at a time
        print(f"\n  Galaxy: {dirname}")
        dirpath = os.path.join(base_folder, dirname)
        file_path = os.path.join(dirpath, filename)
        if os.path.isdir(dirpath) and os.path.isfile(file_path):
            with open(file_path, 'r') as f:
                lines = f.readlines()
                if not header_written and lines:
                    # Write header from the first file
                    all_lines.append(lines[0])
                    header_written = True
                # Add data lines (skip header)
                all_lines.extend(lines[1:])
                rMGMS_galaxy_dic[dirname] = lines[1:]
                # --- Non-parametric comparison to Lin+2020 ---
                xi = np.genfromtxt(lines[1:], usecols=0, dtype=float) # Sigma_star
                xi_mean = np.mean(xi)
                yi = np.genfromtxt(lines[1:], usecols=1, dtype=float) # Sigma_H2
                yi_mean = np.mean(yi)
                delta_xi_mean = np.log10(xi_mean) - np.log10(lin.inverse(yi_mean)) # dex
                delta_yi_mean = np.log10(yi_mean) - np.log10(lin.direct(xi_mean))
                print(f"        Σ_*: {delta_xi_mean:.3f} dex") # physical units (Sol_mass kpc^-2) are omited
                print(f"      Σ_H2 : {delta_yi_mean:.3f} dex")
                galaxies_offset_star_dex[dirname] = delta_xi_mean
                galaxies_offset_h2_dex[dirname] = delta_yi_mean

    # --- Write all collected lines to the output file ---
    with open(output_file, 'w') as f:
        f.writelines(all_lines)
    print(f"\n- Compiled all rMGMSOutputs.txt to: {output_file}")

    # --- Fit rMGMS for all galaxies combined ---
    # all_lines[0] is the header, all_lines[1:] are the data rows
    data = np.genfromtxt(all_lines[1:], dtype=float)
    all_x = data[:, 0]
    all_y = data[:, 1]
    # Plot all the galaxies together with the same color
    gradient, intercept, gradient_err, intercept_err = plotter.plot_rMGMS(a, b, sig, ro, newpath, all_x, all_y, 0) # all_x and all_y are lists, pivot set to 0
    all_x_mean = np.mean(all_x)
    # Plot with a pivot y = m(x-x_pivot) + c is more stable than y = mx + c
    gradient_pivot, intercept_pivot, gradient_pivot_err, intercept_pivot_err = plotter.plot_rMGMS(a, b, sig, ro, newpath, all_x, all_y, all_x_mean) # same as above but using a x_pivot = all_x_mean
    # Plot all the galaxies together using a different color for each (use a list for each one)
    plotter.plot_rMGMS(a, b, sig, ro , newpath, rMGMS_galaxy_dic, 0, 0) # rMGMS_galaxy_dic is a dictionary 0 to funtion as false

    # --- Non-parametric comparison to Lin+2020 - calculate offsets ---
    # Over all galaxies together
    print("\n- Σ_* offset in relation to Lin+2020 rMGMS over all galaxies...")
    all_y_mean = np.mean(all_y)
    delta_xi_total_mean = np.log10(all_x_mean) - np.log10(lin.inverse(all_y_mean))
    delta_yi_total_mean = np.log10(all_y_mean) - np.log10(lin.direct(all_x_mean))
    print(f"            : {delta_xi_total_mean:.3f} dex (Sol_mass kpc^-2)")
    print(f"       Σ_H2 : {delta_yi_total_mean:.3f} dex (Sol_mass kpc^-2)")

    # --- Properties to compare with the offsets ---
    # Sort by galaxy name (dictionary key) so it matches galaxies_offset_star_dex
    stellar_vel_dispersion = [galaxies[gal_name].get('log_stel_vel_disp') 
                            for gal_name in sorted(galaxies.keys()) 
                            if 'log_stel_vel_disp' in galaxies[gal_name]]
    log_SFR = [galaxies[gal_name].get('log_SFR') 
            for gal_name in sorted(galaxies.keys()) 
            if 'log_SFR' in galaxies[gal_name]]
    eff_radius = [galaxies[gal_name].get('log_eff_radius') 
                  for gal_name in sorted(galaxies.keys()) 
                  if 'log_eff_radius' in galaxies[gal_name]]
    log_mu_star = [galaxies[gal_name].get('log_mu_star') 
                   for gal_name in sorted(galaxies.keys()) 
                   if 'log_mu_star' in galaxies[gal_name]]
    properties = {
        "log_stellar_vel_dispersion": [stellar_vel_dispersion, "km/s"],
        "log_SFR": [log_SFR, "Sol_mass year^-1"],
        "log_eff_radius": [eff_radius, "arcsec"],
        "log_mu_star": [log_mu_star, "Sol_mass kpc^2"]
    }

    # --- Statistical tests of correlation between offsets and properties - writting StatisticalTestsOfCorrelation.txt as output ---
    with open(f'{newpath}/StatisticalTestsOfCorrelation.txt', 'w') as f:
        f.write(f"# Correlation coefficients and p-values\n")
    print("\n- Statistical tests of correlation...")
    for key in properties:
        property_values = properties[key][0]
        property_units = properties[key][1]
        galaxies_offset_star_dex_list = list(galaxies_offset_star_dex.values())
        galaxies_offset_h2_dex_list = list(galaxies_offset_h2_dex.values())
        if len(property_values) != len(galaxies_offset_star_dex_list):
            print(f"Warning: Length mismatch between galaxies_offset_star_dex_list and {key}. Skipping this property.")
            continue
        if len(property_values) != len(galaxies_offset_h2_dex_list):
            print(f"Warning: Length mismatch between galaxies_offset_h2_dex_list and {key}. Skipping this property.")
            continue
        # Statistical tests for Σ_star offset
        print(f"\nCorrelations between Σ_star offset with Lin (log) and {key}:")
        pearson_r, pearson_p, spearman_rho, spearman_p, kendall_tau, kendall_p = Statistics.statistical_tests(galaxies_offset_star_dex_list, property_values)
        plotter.plot_parameters(galaxies_offset_star_dex_list, r'$\Sigma_{*}$ Mean Offset from Lin (2020)', "dex", property_values, key, property_units, newpath)
        # Write to file
        with open(f'{newpath}/StatisticalTestsOfCorrelation.txt', 'a') as f:
            f.write(f"\n- Correlations between Σ_star offset with Lin (log) and {key}:\n")
            f.write(f"Pearson: r = {pearson_r:.3f}, p = {pearson_p:.3f}\n")
            f.write(f"Spearman: ρ = {spearman_rho:.3f}, p = {spearman_p:.3f}\n")
            f.write(f"Kendall: τ = {kendall_tau:.3f}, p = {kendall_p:.3f}\n")
        # Statistical tests for Σ_H2 offset
        print(f"\nCorrelations between Σ_H2 offset with Lin (log) and {key}:")
        pearson_r, pearson_p, spearman_rho, spearman_p, kendall_tau, kendall_p = Statistics.statistical_tests(galaxies_offset_h2_dex_list, property_values)
        plotter.plot_parameters(galaxies_offset_h2_dex_list, r'$\Sigma_{H2}$ Mean Offset from Lin (2020)', "dex", property_values, key, property_units, newpath)
        # Write to file
        with open(f'{newpath}/StatisticalTestsOfCorrelation.txt', 'a') as f:
            f.write(f"\n- Correlations between Σ_H2 offset with Lin (log) and {key}:\n")
            f.write(f"Pearson: r = {pearson_r:.3f}, p = {pearson_p:.3f}\n")
            f.write(f"Spearman: ρ = {spearman_rho:.3f}, p = {spearman_p:.3f}\n")
            f.write(f"Kendall: τ = {kendall_tau:.3f}, p = {kendall_p:.3f}\n")
    print(f"\n- Saved correlations to {newpath}/StatisticalTestsOfCorrelation.txt")

    return gradient, intercept, gradient_err, intercept_err, gradient_pivot, intercept_pivot, gradient_pivot_err, intercept_pivot_err


if __name__ == "__main__":
    cli()
