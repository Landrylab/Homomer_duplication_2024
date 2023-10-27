
# # Homodimer identifier
# 
# This script goes through a folder of PDB files and saves suitable homodimers for the simulations. Criteria used to select homodimers:
# - Author assigned biological unit should be "DIMERIC"
# - Experimental data must be crystal structure
# - The two subunits in the file should be identical:
#     - Take the sequences of the subunits from the SEQRES records and check if they are identical (missing residues still appear here)
#     - If the biological unit is a dimer but the file only contains one subunit, check if the second one is derived from transformations (REMARK 350).
#     - Check missing residues section. If the two subunits have different missing residues, discard structure.
#     
#     
# The output will be the folder with the filtered structures and a table mentioning the biological assembly that should be generated for each of them.
#     
# These parameters could also be modified to identify other kinds of complexes
# 

# Load libraries
import glob
import os
import re
from collections import OrderedDict
import shutil
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio import SeqIO

# Define a function to parse pdb lines more easily
def parse_pdb_line(pdb_line):
    '''This function will receive a line from a PDB file and parse it as a list. It will do so based on the
    PDB format explanation from this site:

    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html.
    '''
    atom = pdb_line[0:4].strip(' ')
    atom_num = pdb_line[6:11].strip(' ')
    atom_name = pdb_line[12:16].strip(' ')
    resname = pdb_line[17:20].strip(' ')
    chain = pdb_line[21]
    res_num = pdb_line[22:26].strip(' ')
    x = pdb_line[30:38].strip(' ')
    y = pdb_line[38:46].strip(' ')
    z = pdb_line[46:54].strip(' ')

    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z]


# ## Run on the whole PDB to find all possible homodimers

# Update the variables
structure_type = 'DIMERIC'
pdb_folder = '../../Data/PDB/All_structures'
out_folder = '../../Data/PDB/Homodimers'

aa_three2one = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
    'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G', 
    'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
    'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
    'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'
}

# List all the pdb files in the folder
pdb_files = glob.glob(os.path.join(pdb_folder, '*pdb'))

# Run the loop
# Create the output dictionary if it does not exist
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

# A dictionary for storing the correct biological assembly for each file
# One dictionary for all files, not reinitialized
table_dict = OrderedDict()

# Loop through the PDB files
for pdb_file in pdb_files:

    handle = open(pdb_file, 'r')

    # This boolean will indicate if we have already seen the author determined biological assembly
    # Updated for every file
    bool_author = False

    # A boolean confirming if this is a homodimer or not
    # Updated for every file
    bool_homodimer = False
    
    # Booleans to make sure it is a crystal structure and an adequate length
    bool_length = False
    bool_xtal = False
    
    # A dictionary for the sequences
    # Updated for every file
    seq_dict = {}

    # A list for the number of biomolecular transformations that should be applied
    # Updated for every file
    biomt_list = []

    # A dictionary for missing residues
    # Updated for every file
    misres_dict = OrderedDict()

    for pdb_line in handle:

        # Check the type of experimental data
        if pdb_line.startswith('EXPDTA'):
            # Two ifs make sure I am looking at the correct line when deciding to skip
            if pdb_line.startswith('EXPDTA    X-RAY DIFFRACTION'):
                bool_xtal = True
            else:
                # If this is not a crystal structure, skip this structure
                break
        
        ## Get the author determined biological assembly and transformations to generate other subunits (REMARK 350)
        if pdb_line.startswith('REMARK 350') and not bool_author:

            # Save the biological assembly we are looking at
            temp_bio = re.search(pattern = 'BIOMOLECULE:\s+(\d+)', string = pdb_line)

            if temp_bio:
                bio_assembly = temp_bio.group(1)

            temp_struc = re.search(pattern = 'AUTHOR DETERMINED BIOLOGICAL UNIT:\s+([A-Za-z]+)', string = pdb_line)

            if temp_struc:
                bio_struc_type = temp_struc.group(1)
                bool_author = True

        ## Check for transformations to generate other subunits once we have the author determined
        ## biological assembly (REMARK 350)
        if pdb_line.startswith('REMARK 350') and bool_author:
            # Look for rows that contain BIOMT
            if re.search(pattern = 'BIOMT\d+', string = pdb_line):
                biomt_line = re.split(pattern = '\s+', string = pdb_line)
                if not int(biomt_line[3]) in biomt_list:
                    biomt_list.append(int(biomt_line[3]))


        ## Get the sequence for the subunits from the SEQRES
        if pdb_line.startswith('SEQRES'):
            seqres_line = re.split(pattern = '\s+', string = pdb_line)
            subunit = seqres_line[2]

            # Make sure this subunit is already in the dictionary
            if seq_dict.get(subunit, -1) == -1:
                seq_dict[subunit] = ''

            for residue in seqres_line[4:]:
                if residue != '':
                    # Residue types have to be included in the list of amino acid types
                    # This will remove structures containing nucleic acids and structures with modified
                    # residues
                    if aa_three2one.get(residue , -1) == -1:
                        bool_homodimer = False
                        continue
                    else:
                        seq_dict[subunit] = seq_dict[subunit] + aa_three2one[residue]

        ## Check the missing residues (REMARK 465)
        if pdb_line.startswith('REMARK 465'):
            # Split by whitespace
            misres_line = re.split(pattern = '\s+', string = pdb_line)

            # Check length equals 6 (considering the empty string at the last position)
            if len(misres_line) == 6:
                # Position 3 (0-based) is the chain
                subunit_misres = misres_line[3]

                if not subunit_misres in misres_dict.keys():
                    misres_dict[subunit_misres] = []

                # Save this residue to the list of missing residues for this chain
                missing_residue = misres_line[2] + misres_line[4]
                misres_dict[subunit_misres].append(missing_residue)

    
    ## Use all the data to determine if it is a homodimer
    # If there is no author determined biological assembly, move on to the next structure
    if not bool_author:
        continue
    
    # Check if author determined assembly is a dimer
    if bio_struc_type == structure_type:
        # Check if there is only one subunit in the file from which the second one will be generated
        if len(seq_dict.keys()) == 1 and max(biomt_list) == 2:
            bool_homodimer = True

        # Check if the two subunits in the file are identical
        elif len(seq_dict.keys()) == 2 and list(seq_dict.values())[0] == list(seq_dict.values())[1]:
            bool_homodimer = True

            # Check if missing residues are the same for the two chains
            if len(misres_dict.values()) == 2:
                if list(misres_dict.values())[0] != list(misres_dict.values())[1]:
                    bool_homodimer = False
            # Missing residues are not the same for the two chains if there are only missing residues on one
            # chain. 
            elif len(misres_dict.values()) == 1:
                bool_homodimer = False


    handle.close()

    # Parse the pdb file
    chain = {record.id: record.seq for record in SeqIO.parse(pdb_file, 'pdb-seqres')}
    
    # Handle structures like 1hya that only have HETATM, so no chains
    if len(list(chain.keys())) > 0:
        query_chain = chain[list(chain.keys())[0]]
    else:
        # Files with no chains can be discarded
        continue

    # Get the length of the chain
    sequence_length = len(query_chain)
    
    if sequence_length > 60 and sequence_length < 450:
        bool_length = True
    
    if bool_homodimer and bool_xtal and bool_length:
        # Copy file to output directory
        shutil.copy(pdb_file, out_folder)

        # Add to a dictionary with the PDB files and the correct bio assembly for each
        file_basename = os.path.splitext(os.path.basename(pdb_file))[0]
        table_dict[file_basename] = [file_basename, bio_assembly, sequence_length]

## After the loop
# Convert the dictionary to a data frame and save as a table
assembly_dataframe = pd.DataFrame.from_dict(table_dict, columns = ['Structure', 'Bio_assembly', 'Length'], orient = 'index')
assembly_dataframe.to_csv(path_or_buf = os.path.join(out_folder, 'bio_assembly_table.txt'), index = False, 
                          sep = '\t')

pdb_file


# ## Now use the file produced in the above loop to extract the sequences from the FASTA file with the whole PDB
# 

# Load a list with all the files from the folder
pdb_list = glob.glob('../../Data/PDB/Homodimers/*.pdb')

# Extract the PDB IDs
pdb_list = [entry.split('/')[-1][0:4] for entry in pdb_list]

print(pdb_list[0:3])
print(len(pdb_list))

pdb_dict = {}

# Convert the list to a dictionary
for entry in pdb_list:
    pdb_dict[entry] = 1

# Open the fasta file with all the sequences
whole_pdb = SeqIO.parse('../../Data/PDB/pdb_seqres_2023-03-29.txt', 'fasta')

records = []

for entry in whole_pdb:
    
    # Split the description to make sure we are looking at proteins
    metadata = entry.description.split(' ')
    seq_moltype = metadata[1]
    
    seq_id = entry.id[0:4]
    
    # Check if this entry is in the dictionary
    if pdb_dict.get(seq_id, -1) != -1:
        # Check if the entry has yet to be added and make sure the moltype is a protein
        if pdb_dict[seq_id] == 1 and seq_moltype == 'mol:protein':
            # Add the entry to a list of records to be written to a new file
            records.append(entry)
            
            # Mark this entry in the dictionary so that we don't add it again (since they
            # are homodimers I only need one copy of the sequence)
            pdb_dict[seq_id] = 0
            

len(records)

list_records = []
for entry in records:
    list_records.append(entry.id[0:4])
    
for entry in pdb_list:
    if not entry in list_records:
        print(entry)


# Records 5dqx and 7lj6 are obsolete entries, so they were discarded.
# Write the records to a new file
SeqIO.write(records, '../../Data/PDB/pdb_seqres_2023-03-29_filtered.txt', 'fasta')


# # CD-Hit analysis
# 
# This section will use CD-Hit to cluster the 104 final structures by sequence identity.
## Use the following command  to generate the CD-Hit clusters


## Load CD-HIT
module load cd-hit/4.8.1

## Run CD-Hit with all the sequences from the PDB as of March 29th, 2021
# -i input file
# -o output file prefix (fasta representative sequences, text file of clusters)
# -n word size (3 is recommended in the user guide for thresholds of 0.5 - 0.6)
# -c sequence identity threshold (0.5)
# -M max available memory in Mb 
# -G global alignment (1) or local alignment (0)
# -s length difference cutoff (sequences have to be at least this fraction of the length of the representative to be in the cluster)
# -T number of threads

cd-hit -i ../Data/pdb_seqres_filtered_new.txt -o Results_c0.4_new/pdb_seqres_filtered_cdhit_40pct_new -n 2 -c 0.4 -M 8000 -G 1 -s 0.8 -T 5

# Load libraries
import numpy as np
import pandas as pd
import os
import sys
import re

## Read the list of PDB structures
list_structures = pd.read_csv('../../Data/final_104_structures.txt', names = ['Complex'])
list_structures


## Load the CD-Hit clustering results (40% sequence identity)
infile = '../../Data/CD-Hit/Results/Results_c0.4/pdb_seqres_filtered_cdhit_40pct.clstr'

pdb2cluster_40pct = {}
cluster = 0
with open(infile, 'r') as handle:
    for line in handle:

        # Check current cluster
        if line.startswith('>'):
            cluster = int(re.search(string = line, pattern = '>Cluster (\d+)').group(1))
        else:
            # Extract all the protein sequences
            pdb_id = re.search(string = line, pattern = '>(.*)_').group(1)

            # Add each to the dictionary
            pdb2cluster_40pct[pdb_id] = cluster
            
list_structures['cluster_40pct'] = [pdb2cluster_40pct.get(pdb, 'NA') for pdb in list_structures['Complex']]

list_structures

list_structures['cluster_40pct'].value_counts()

