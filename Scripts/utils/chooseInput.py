from InquirerPy import inquirer
import os
import shutil
from .sideQuests import getValues

# --- General user input handling class ---
class userChoice:
    def __init__(self, options):
        self.options = options

    def galaxyChoice(self):
        choice = inquirer.select(
            message="\n> Select a galaxy:",
            choices=self.options,
        ).execute()
        print(f"You selected: {choice}")

        return choice

    def runChoice(self):
        choice = inquirer.select(
            message="n> Select the name of the run (folder name):",
            choices=self.options,
        ).execute()
        print(f"You selected: {choice}")

        return choice

    # This one is exclusive to runs only, excludes runsets   
    def runExclusiveChoice(self):
        self.options = [s for s in self.options if not s.startswith("RUNSET")]
        choice = inquirer.select(
            message="n> Select the name of the run (folder name):",
            choices=self.options,
        ).execute()
        print(f"You selected: {choice}")

        return choice
    
    def runsToCompare(self):
        i = 0
        run_name = 0
        chosen_runs = []
        while run_name != "Continue":
            i += 1
            run_name = inquirer.select(
                message="\n> Select the runs to use (folder name):",
                choices=self.options,
            ).execute()
            del self.options[self.options.index(run_name)]
            chosen_runs.append(run_name)
            if (i==1) and (not run_name.startswith("RUNSET")):
                # Filter the list to keep only those that do NOT start with RUNSET
                self.options = [s for s in self.options if not s.startswith("RUNSET")]
            if i == 2: 
                self.options.append("Continue")
            if (i==1) and (run_name.startswith("RUNSET")):
                self.options = [s for s in self.options if not s.startswith("RUNSET")]
                chosen_set = f"./Runs/{run_name}"
                chosen_runs = [f"{run_name}/" + name for name in sorted(os.listdir(chosen_set)) if os.path.isdir(os.path.join(chosen_set, name))]
                run_name = "Continue"
            elif (i > 2) and (run_name=="Continue"):
                del chosen_runs[chosen_runs.index(run_name)]
        print("\nChosen runs:")
        [print(f'  {x}') for x in chosen_runs]

        return chosen_runs
    
    def whatToDo(self):
        scriptToRun = inquirer.select(
            message="\n> What to do...",
            choices=self.options,
        ).execute()

        return scriptToRun
    
    def deleteRuns(self):
        self.options.append("Exit")
        run_to_delete = 0
        while run_to_delete != "Exit":
            run_to_delete = inquirer.select(
                message="\n> Select a run to delete:",
                choices=self.options,
            ).execute()
            if run_to_delete!="Exit":
                getValues.get_folder_size('Runs/' + run_to_delete)
                delete_sure = input(f"Are you sure you want to delete {run_to_delete}? (y/n): ").strip().lower()
                if delete_sure == "y":
                    del self.options[self.options.index(run_to_delete)]
                    shutil.rmtree(f"Runs/{run_to_delete}")
                    
# --- Run initial input handling class ---
class initialInput:
    def __init__(self):
        pass

    def readMe():
        readme = input("\nRead README.md? (y/n): ").strip().lower()
        if readme == "y":
            print("")
            with open("README.md", "r") as file:
                content = file.read()
                print(content)
        else: 
            print("- Skipping README.md")

    def speedRun():
        show_figures = input("\nShow figures during this run? (y/n): ").strip().lower()
        if show_figures == "y":
            show_figures = True
            print("- Figures will be shown during this run. This may slow down the process.")
            fast_run = False
        else:
            show_figures = False
            print("- Figures will NOT be shown during this run. This will speed up the process.")
            fast_run = input("\nSkip unnecessary user prompts to fasten the process even further? (y/n): ").strip().lower()
            if fast_run == "y":
                fast_run = True
            else:
                fast_run = False

        return show_figures, fast_run
    
    def matchMaps(lowest_resolution):
        match_maps = input(f"\nMatch all maps to a specific resolution? (y/n/{lowest_resolution[0]}'s(g)): ").strip().lower()
        if match_maps == "g":
            print(f"- All maps will be matched to the lowest resolution: {lowest_resolution[0]} with {lowest_resolution[1]:.2f} pc/pixel\n")
            match_maps = True
        elif match_maps == "y":
            resolution_input = input("Input a resolution (in parsec): ")
            lowest_resolution[1] = float(resolution_input)
            match_maps = True
        else:
            print("- Maps will NOT be matched to the lowest resolution.\n")
            match_maps = False
        
        return match_maps
