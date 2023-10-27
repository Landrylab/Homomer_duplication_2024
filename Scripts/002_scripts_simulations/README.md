Scripts in this folder include:

- 008\_simulation\_workflow.sh: This is the script that manages the simulations. Some input parameters have to be modified inside the script to execute the different regimes of simulations. When running parametric simulations, only the first two arguments (boolean to confirm parametric simulations and results folder) are required. When running simulations with PDB structures, arguments 3 and 4 (path to final matrices of mutational effects, PDB ID) are required as well. This script then uses the following scripts:

	. input\_params\_simulations\_slurm.py: This file contains the main set of parameters to be loaded in the simulations.

	. 008\_simulations\_slurm\_manager.sh: This script submits the job to a SLURM manager and calls the main simulations script.

	. 008\_simulations\_slurm.py: This is the main script for the simulations.
