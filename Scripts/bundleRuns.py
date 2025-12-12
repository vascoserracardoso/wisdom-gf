#!/usr/bin/env python3

# --- Import packages ---
import os
from InquirerPy import inquirer

# --- Import custom modules ---
from .analysis import Plotter
from .utils import userChoice

# --- Command line interface ---
def cli():
    imported_bundle = True
    chosen_runs = userChoice(options = [name for name in sorted(os.listdir('./Runs')) if os.path.isdir(os.path.join('./Runs', name))]).runsToCompare()
    show_figures = input("\nShow figures? (y/n): ").lower() == 'y'
    main(chosen_runs,
         show_figures)

def main(runs, show_figures):

    quantities = {}
    for run in runs:
        path_info = f'./Runs/{run}/keyValues.txt'
        with open(path_info, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split()

                var, value, unit = parts[0], parts[1], parts[2]

                # Convert to float if possible
                try:
                    value = float(value)
                except ValueError:
                    continue  # skip if not numeric

                # Append to dictionary list
                if var not in quantities:
                    quantities[var] = [unit] # The first element of the list are the units
                quantities[var].append(value)

    Plotter(show_figures).dic_plot(quantities, 'Resolution', 'Gradient', 'Gradient_err')
    Plotter(show_figures).dic_plot(quantities, 'Resolution', 'Intercept', 'Intercept_err')
    Plotter(show_figures).dic_plot(quantities, 'Resolution', 'Gradient_pivot', 'Gradient_pivot_err')
    Plotter(show_figures).dic_plot(quantities, 'Resolution', 'Intercept_pivot', 'Intercept_pivot_err')


if __name__ == "__main__":
    cli()

