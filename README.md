Scripts and data for Cisneros et al., (2023) Mutational biases contribute to the neutral complexification of protein interaction networks following gene duplication.

Folders contain the following:
- Data: It contains folders to store the raw PDB files and intermediate files created by the data preparation steps, such as
the CD-Hit sequence clustering, the DSSP secondary structure analysis, and the matrices with the distributions of mutational
effects. The following data should be downloaded from their respective databases:

	. ECOD data: Should be stored in the Data/Structures folder. The ecod_latest.domains file should be downloaded from:
	http://prodata.swmed.edu/ecod/
	

	. PDB data: All the structure files should be stored in the Data/PDB/All\_structures folder. Structures and the pdb_seqres.txt file should be downloaded from:
	https://www.rcsb.org/docs/programmatic-access/file-download-services#pdb-entry-files

A helper repository with large files should be downloaded from Zenodo (https://zenodo.org/records/10048861).
The data from this repository should be uncompressed and organized as follows:
	- all\_mut\_matrices.tsv: Should be uncompressed and moved to Homomer\_duplication\_2023/Data
	- Results_pdup1: Should be uncompressed and moved to Homomer\_duplication\_2023/Results\_simulations
	- Results\_solution\_space: Should be uncompressed and moved to Homomer\_duplication\_2023
	- Structures: Should be uncompressed and moved to Homomer\_duplication\_2023/Data

- Figures: pdf versions of the main and supplementary figures

- Scripts: Contains scripts for data preparation, analysis, and figure generation. Folders and scripts are numbered according
to the order in which they were used in the analysis.

- Results\_simulations: Use this folder to store the results of the simulations. This folder will contain the Results\_pdup1 from Zenodo
once it is decompressed. The two subfolders inside contain the results of the different simulations for both the parametric and single
structure regimes.

- Results\_solution\_space: Use this folder to store the results of our exploration of the solution space of the system
of equations. These can be downloaded from Zenodo and decompressed. 


