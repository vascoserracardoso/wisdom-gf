# WISDOM-GF

WISDOM-GF is a software package developed to investigate the resolved Molecular Gas Main Sequence (rMGMS) of galaxies at parsec-scale resolution. It computes gas and stellar mass surface densities, performs pixel-by-pixel comparisons between the resulting maps, and derives the rMGMS relation together with other relevant physical quantities. Each galaxy is analysed within its own dedicated directory, which contains the required FITS inputs and the output products generated automatically by the program. Once the software is running, the workflow is largely self-explanatory, and the resulting files and figures follow a clear structure that facilitates interpretation.

Requirements and installation:
WISDOM-GF requires a standard scientific Python environment. All necessary dependencies can be installed automatically using the provided Makefile (previous requirements: Python3, pip and Python venv module). Running the command “make” in the project’s root directory creates a virtual environment and installs all required packages. In case the creation of the virtual environment fails, dependencies are installed directly on the system. The full analysis can then be executed using “make run”. To remove the virtual environment, use “make clean”. Alternatively, users who prefer to manage dependencies manually may run “python main.py” from the root directory, provided all necessary packages are installed.

Running the software:
The complete analysis is executed using either “python main.py” or “make run”. The software processes all steps sequentially, assuming that the directory structure remains unchanged. Individual scripts in the Scripts directory may also be executed manually for development or debugging, but they must always be run from the root directory due to relative path dependencies.

Adding a new object:
To analyse a new galaxy, copy the NGCexample folder into the Galaxies directory and rename it to match the target object. Place the moment-0 FITS map inside Galaxies/GalaxyName/FITS/inputs using the filename GalaxyName_mom0.fits, and add the corresponding parameters to Scripts/Dictionary.py. All other required files (including processed FITS maps, figures, and summary outputs) are generated automatically during execution. If the FITS header does not include the rest frequency or other required metadata, these values must be provided manually in Dictionary.py.

Data availability:
This repository does not distribute any observational data. The FITS files required for the analysis remain the property of their respective providers and must be obtained independently by the user. The software will not run without these data. The MIT License included in this repository applies only to the source code and documentation, not to any external data products.

Outputs:
For each analysed galaxy, the software produces processed FITS files, diagnostic figures, and a summary of physical quantities. These outputs follow consistent naming conventions and are stored within the corresponding galaxy’s directory. The organisation of the output files makes it straightforward to follow the analysis steps and interpret the results.

Contact and citation:
If you use WISDOM-GF in academic work, please cite this repository or acknowledge its use appropriately. For questions, suggestions, or issues related to the software, please contact the author or submit an issue through the repository.
