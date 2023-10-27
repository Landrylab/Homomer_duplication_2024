Scripts in this folder include:

- 000\_homodimer\_identifier.py: This is the script that surveys the PDB for suitable homodimers for the simulations. It also checks the CD-Hit clustering to make sure that the selected structures are not homologous.

- 001\_generate\_bio\_assembly.py: This script generates the author-designated biological assembly for each PDB file.

- 002\_insertion\_residues\_alternative\_coordinates.py: This script renumbers insertion residues and keeps coordinates with the highest occupancy for each residue. This is useful when simulating mutations with FoldX to help it recognize each residue correctly.

- 003\_foldx\_repair\_slurm.sh and 003\_foldx\_repair\_slurm\_multiple.sh: Scripts used to run FoldX\_Repair 10 times on a single structure to minimize free energy by solving bad contacts.

- 004\_call\_interfaces.py and call\_interfaces\_helper.py: The script that identifies the interfaces of the homodimers.

- 005\_mutatex\_slurm.sh: Script to call MutateX and generate the distributions of mutational effects on folding energy and binding affinity.

- 006\_gather\_mutations.sh: Script to gather the MutateX generated data.

- 007\_mutation\_parameters\_Manitou.R: Script to organize the gathered MutateX data into a matrix of mutational effects.

- 1aof\_sasa.pdb: A reference PDB file to get van der Waals radii for different atom types.

