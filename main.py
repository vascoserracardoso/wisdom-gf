#!/usr/bin/env python3

# --- Import packages ---
import os
os.system('cls' if os.name == 'nt' else 'clear')
print("Loading packages, scripts and modules...")
from datetime import datetime
import numpy as np
import shutil
from astropy.io import fits
import subprocess
import sys

# --- Import scripts ---
import Scripts.mom0_to_SigmaH2 as mom0
import Scripts.mgeParams_to_SigmaStar as mge
import Scripts.compareSigmas as comp
from Scripts.Dictionary import galaxies
import Scripts.addrMGMS as addrMGMS
import Scripts.bundleRuns as bundle

# --- Import modules ---
from Scripts.utils import initialInput
from Scripts.utils import getValues
from Scripts.utils import userChoice

# --- Define analysis function ---
def analysis(galaxy, show_figures, match_maps, lowest_resolution, fast_run):

    # --- Initialize PhysicalQuantitiesOutputs.txt file ---
    with open(f'./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt', 'w') as f:
        f.write("# Physical quantities calculated for the galaxy " + galaxy + "\n")

    # --- Run mom0_to_SigmaH2.py and mgeParams_to_SigmaStar.py scripts ---
    import_scripts = True
    print("\n> Running mom0_to_SigmaH2.py...\n")
    mom0.main(galaxy, import_scripts, show_figures, match_maps, lowest_resolution)
    if fast_run == False:
        input("\nPress Enter to continue...")
    print("\n> Running mgeParams_to_SigmaStar.py...\n")
    mge.main(galaxy, import_scripts, show_figures, match_maps, lowest_resolution)
    if fast_run == False:
        input("\nPress Enter to continue...")

    # --- Total mass computation, gas fraction calculation and writing to file ---
    print("\nTotal masses and gas fraction...")
    physicalQuantities = {}
    with open(f'./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt', 'r') as f:
        next(f) # Ignore first line
        for line in f:
            name, value, unit = line.split()
            physicalQuantities[name] = float(value)
    M_H2 = physicalQuantities['M_H2']
    M_star = physicalQuantities['M_star']
    print(f"- M_H2 = {M_H2:.6e} {unit}")
    print(f"- M_star = {M_star:.6e} {unit}")
    gas_fraction = M_H2 / (M_star + M_H2)
    print(f"- Gas fraction: {gas_fraction:.2%} of total mass")
    M_H2_reference = galaxies[galaxy]['H2_mass']
    M_star_reference = galaxies[galaxy]['stellar_mass']
    print(f"- M_H2(from reference) = {M_H2_reference:.6e} {unit}")
    print(f"- M_star(from reference) = {M_star_reference:.6e} {unit}")
    M_H2_er = np.abs(M_H2 - M_H2_reference) / M_H2_reference
    M_star_er = np.abs(M_star - M_star_reference) / M_star_reference
    print(f"- Error in M_H2 compared to reference: {M_H2_er:.2%}")
    print(f"- Error in M_star compared to reference: {M_star_er:.2%}")
    with open(f'./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt', 'a') as f:
        f.write(f"Gas_fraction        {gas_fraction:.2%} no_unit\n")
        f.write(f"M_H2(Reference)     {M_H2_reference:.6e} {unit}\n")
        f.write(f"M_star(Reference)   {M_star_reference:.6e} {unit}\n")
        f.write(f"Error(M_H2)         {M_H2_er:.2%} no_unit\n")
        f.write(f"Error(M_star)       {M_star_er:.2%} no_unit\n")
    print(f"- Saved total masses from reference to ./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt")

    if fast_run == False:
        input("\nPress Enter to continue...")

    # --- Run compareSigmas.py script ---
    print("\n> Running compareSigmas.py...\n")
    comp.main(galaxy, import_scripts, show_figures, match_maps, lowest_resolution)


def main(scale):
    # --- Create folder to save outputs of each run ---
    run_folder_name = datetime.now().strftime("%Y%m%d_%H%M%S")
    path_run_folder = './Runs/' + run_folder_name
    os.makedirs(path_run_folder, exist_ok=True)
    print('\nRun outputs will be saved to ', path_run_folder)

    input("\nPress Enter to start...")

    # --- Get user input regarding run options ---
    show_figures, fast_run = initialInput.speedRun()

    print("\n> Checking galaxies dictionary...\n")

    # --- Check the maximum number of keys in the galaxies dictionary ---
    max_keys = max(len(subdict) for subdict in galaxies.values())
    for g, d in galaxies.items():
        print(f"{g}: {len(d)} keys")
    print(f"Max keys = {max_keys}. Only galaxies with {max_keys} keys will be analysed.\n")

    # --- Get lowest resolution among galaxies with complete set of keys ---
    lowest_resolution = getValues.lowestResolution(galaxies, max_keys)
    print(f"Lowest resolution in the set of galaxies: {lowest_resolution[0]} with {lowest_resolution[1]:.2f} pc/pixel")
    if scale !=0:
        match_maps = True
        lowest_resolution[1] = float(scale)
    else:
        match_maps = initialInput.matchMaps(lowest_resolution)
    print(f"\nResolution in use: {lowest_resolution[1]}")

    # --- Run analysis for all galaxies with the maximum number of keys to discard incomplete entries ---
    for galaxy_name, subdict in galaxies.items():
        if len(subdict) == max_keys:
            print(f'\n\n######################################')
            print(f'\n\n##### Analyses of galaxy {galaxy_name} #####')
            print(f'\n\n######################################')
            folder_clear = f"./Galaxies/{galaxy_name}/FITS/outputs" 
            print(f'\n> Deleting all FITS files from {folder_clear} <')
            for filename in os.listdir(folder_clear):
                file_path = os.path.join(folder_clear, filename)
                try:
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                        print(f"Deleted file: {file_path}")
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")
            analysis(galaxy_name, show_figures, match_maps, lowest_resolution, fast_run)
            print(f'\n\n#######################################')
            print(f'\n\n##### End of analyses for {galaxy_name} #####')
            print(f'\n\n#######################################')

            if fast_run == False:
                input('\nPress Enter to continue...')


    # --- Copy the galaxy folder to the run folder ---
    src_folder = f'./Galaxies'
    dst_folder = f'{path_run_folder}/Galaxies'
    shutil.copytree(src_folder, dst_folder, dirs_exist_ok=True)
    print(f"- Copied {src_folder} to {dst_folder}")

    # --- Run addrMGMS.py script as imported script ---
    imported_script = True
    print("\n> Running addrMGMS.py...\n")
    gradient, intercept, gradient_err, intercept_err, gradient_pivot, intercept_pivot, gradient_pivot_err, intercept_pivot_err = addrMGMS.main(imported_script,run_folder_name,show_figures)

    if fast_run == False:
        input("\nPress Enter to continue...")

    # --- Save general information to keyValues.txt file ---
    with open(f'{path_run_folder}/keyValues.txt', 'a') as f:
        if match_maps:
            f.write("# Important inputs and outputs of this run\n")
            f.write(f"Resolution        {lowest_resolution[1]:.2f} pc\n")
            f.write(f"Gradient         {gradient:.5f} no_unit\n")
            f.write(f"Gradient_err         {gradient_err:.5f} no_unit\n")
            f.write(f"Intercept        {intercept:.5f} dex\n")
            f.write(f"Intercept_err        {intercept_err:.5f} dex\n")
            f.write(f"Gradient_pivot         {gradient_pivot:.5f} no_unit\n")
            f.write(f"Gradient_pivot_err         {gradient_pivot_err:.5f} no_unit\n")
            f.write(f"Intercept_pivot        {intercept_pivot:.5f} dex\n")
            f.write(f"Intercept_pivot_err        {intercept_pivot_err:.5f} dex\n")
        else:
            f.write(f"Resolution:        None no_unit\n")

    # --- Rename run folder to include the scale used get the its respective size ---
    os.rename(path_run_folder, path_run_folder + f"_SCALE_{lowest_resolution[1]:.2f}")
    path_run_folder = path_run_folder + f"_SCALE_{lowest_resolution[1]:.2f}"
    print(f"\nRun folder name updated to {path_run_folder}")
    print(f"- Saved general information to .{path_run_folder}/keyValues.txt'")
    getValues.get_folder_size(path_run_folder)

    input("\nPress Enter to exit run...")

    return path_run_folder


if __name__ == "__main__":

    # --- Show readMe information ---
    initialInput.readMe()

    # --- Main program loop ---
    program_on = True
    while program_on:

        # --- Show program banner ---
        print(r"""

         _       ___________ ____  ____  __  ___            ____________
        | |     / /  _/ ___// __ \/ __ \/  |/  /           / ____/ ____/
        | | /| / // / \__ \/ / / / / / / /|_/ /  ______   / / __/ /_    
        | |/ |/ // / ___/ / /_/ / /_/ / /  / /  /_____/  / /_/ / __/    
        |__/|__/___//____/_____/\____/_/  /_/            \____/_/       
                                                                

        """)

        # --- Show size of WISDOM-GF folder ---
        getValues.get_folder_size('../wisdom-gf')

        # --- Check if running inside conda or virtual environment ---
        conda_environment = False
        virtual_environment = False
        if "CONDA_PREFIX" in os.environ: # This checks if any conda environment is on (even base)
            print(f"Running inside a Conda environment!")
            conda_environment = True
        exe = sys.executable
        venv_path = os.path.join(os.getcwd(), ".venv")
        if exe.startswith(venv_path): # This one checks if program was ran using the virtual environment (make run)
            print("Running inside project virtual environment (.venv)")
            virtual_environment = True
        print(f'Python interpreter: {sys.executable}')

        # --- Show options that the program offers and get user choice ---
        scripts_folder = './Scripts'
        options = ['(1) Calculate map of GAS mass surface density (Σ_H2) of an object.', \
                '(2) Calculate map of STELLAR mass surface density (Σ_*) of an object.', \
                '(3) Compare maps of Σ_H2 and Σ_* of an object.', \
                '(4) Construct rMGMS relationship from a previous SINGLE RUN.', \
                '(5) Perform a SINGLE RUN. Run 1, 2, 3 for all available objects and then run 4.', \
                '(6) Perform multiple runs (RUNSET) at different scales.', \
                '(7) Compare set of runs (RUNSET) or single runs at different scales.', \
                '(8) Delete runs.', \
                '(9) Exit.']
        to_do = userChoice(options).whatToDo()
        if to_do == options[0]:
            if virtual_environment:
                subprocess.run([".venv/bin/python", "-m", "Scripts.mom0_to_SigmaH2"])
            else:
                subprocess.run(["python", "-m", "Scripts.mom0_to_SigmaH2"])
        if to_do == options[1]:
            if virtual_environment:
                subprocess.run([".venv/bin/python", "-m", "Scripts.mgeParams_to_SigmaStar"])
            else:
                subprocess.run(["python", "-m", "Scripts.mgeParams_to_SigmaStar"])
        if to_do == options[2]:
            if virtual_environment:
                subprocess.run([".venv/bin/python", "-m", "Scripts.compareSigmas"])      
            else:  
                subprocess.run(["python", "-m", "Scripts.compareSigmas"])
        if to_do == options[3]:
            if virtual_environment:
                subprocess.run([".venv/bin/python", "-m", "Scripts.addrMGMS"])      
            else:  
                subprocess.run(["python", "-m", "Scripts.addrMGMS"])
        if to_do == options[4]:
            scale = 0
            print('\n\n##########################################################################')
            print("\n> main.py running as script\n")
            main(scale)
        # Perform multiple runs at different scales and move them to a runset folder in the end
        if to_do == options[5]:
            scales_input = input('\nEnter scales in parsec (e.g. 20,30,40,... or range:10,100,5): ')
            if scales_input.startswith('range:'):
                range = [float(x) for x in scales_input[6:].split(",")]
                small = range[0]
                big = range[1]
                interval = range[2]
                scales = []
                while small <= big:
                    scales.append(small)
                    small += interval
            else:
                scales = [float(x) for x in scales_input.split(",")]
            folders_to_move = []
            for scale in scales:
                print('\n\n##########################################################################')
                print("\n> main.py running as script\n")
                runName = main(scale)
                folders_to_move.append(runName)
            new_folder_path = "./Runs/RUNSET_" + datetime.now().strftime("%Y%m%d_%H%M%S")
            os.makedirs(new_folder_path, exist_ok=True)
            print(f"\nCreated folder: {new_folder_path}")
            for folder in folders_to_move:
                src = folder
                dst_folder = folder[7:]
                dst = os.path.join(new_folder_path, dst_folder)
                if os.path.exists(src):
                    shutil.move(src, dst)
                    print(f"Moved '{folder}' into '{new_folder_path}'")
                else:
                    print(f"Warning: '{folder}' not found in ./Runs")
            getValues.get_folder_size(new_folder_path)
        if to_do == options[6]: 
            if virtual_environment:
                subprocess.run([".venv/bin/python", "-m", "Scripts.bundleRuns"])      
            else:  
                subprocess.run(["python", "-m", "Scripts.bundleRuns"])
        if to_do == options[7]: 
            userChoice(options = [name for name in sorted(os.listdir('./Runs')) if os.path.isdir(os.path.join('./Runs', name))]).deleteRuns()
        if to_do == options[8]:
            print("All prints will be deleted from the terminal...")
            make_sure = input("Are you sure you want to exit ? (y/n)").strip().lower()
            if make_sure == 'y':
                program_on = False
                os.system('cls' if os.name == 'nt' else 'clear')
            else:
                continue
