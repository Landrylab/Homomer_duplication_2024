# # Summarize DSSP files 
# 
# This script will read a series of dssp output files and generate a summary for each that I could try to relate to the HET / HM ratio in the simulations. For reference:
# 
# - H = alpha-helix
# - B = beta-bridge residue
# - E = extended strand (in beta ladder)
# - G = 3/10-helix
# - I = 5-helix
# - T = H-bonded turn
# - S = bend

# Load libraries
import csv
import numpy as np
import pandas as pd
import re
import glob
import os

code_dict = {
    'H':'Alpha helix',
    'B':'Beta bridge',
    'E':'Beta ladder', 
    'G':'3/10 helix',
    'I':'5-helix',
    'T':'H-bonded turn',
    'S':'Bend'
}

## Define some functions
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
    res_num = pdb_line[22:27].strip(' ')
    x = pdb_line[30:38].strip(' ')
    y = pdb_line[38:46].strip(' ')
    z = pdb_line[46:54].strip(' ')
    bfactor = pdb_line[60:66].strip(' ')

    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z, bfactor]

# ## Need to adapt this block and make it loop through the list of structures

## Folder with the DSSP results
folder_path = '../../Data/DSSP'

# Get the list of files in that folder
file_list = glob.glob(os.path.join(folder_path, '*dssp'))

df_used_structures = pd.read_csv('../../Data/final_104_structures.txt', sep = '\t', names = ['PDB_ID'])
used_structures = [entry for entry in df_used_structures['PDB_ID']]
used_structures

out_folder = os.path.join(folder_path, 'Clean_tables')

os.makedirs(out_folder, exist_ok = True)

path_interfaces = '../../Data/Structures/004_interfaces/'

summary_dict = {}
summary_dict_interfaces = {}

# Loop through the files
for infile in file_list:
    # Extract the pdb id from the file name
    file_basename = os.path.basename(infile)
    pdb_id = str(file_basename.split('_')[0])
    # print(pdb_id)
    
    # Skip if this structure is not in the final data set
    if not pdb_id in used_structures:
        continue
    
    handle = open(infile, 'r')

    # Start a dictionary to keep track of the values
    sec_struc_annotations = {}

    # A boolean to help us know when to start looking at the data
    bool_data = False
    
    # Loop through the lines
    for line in handle:

        # Looking for the start of the data
        if not bool_data:
            if line.startswith('  #'):
                # The hashtag marks the header, so we know that we should start reading data from the 
                # following line
                bool_data = True
                continue
            else:
                continue
        
        else:
            # Start looking at the data from here

            # Split by one or more instances of whitespace
            split_line = re.split(pattern = '\s+', string = line)

            # Columns 2 and 5 (0-based) contain the position and the secondary structure annotation
            # However, positions with no annotation in column 5 will skip that column, so we need to
            # make sure the value in that column is a secondary structure annotation

            # If column 2 is "!*", it means we have finished with the first subunit. I will stop there since 
            # the structure is a symmetric homomer
            if split_line[2] == '!*' or split_line[2] == '!':
                break
            elif split_line[2] == '!':
                continue
            else:
                position = int(split_line[2])
                sec_struc = split_line[5]
            
            # Save solvent accessibility
            solv_acc = int(line[34:39].strip())

            # If this is a secondary structure annotation
            if sec_struc in ('H', 'B', 'E', 'G', 'I', 'T', 'S'):
                sec_struc_annotations[position] = [position, code_dict[sec_struc], solv_acc]
            else:
                sec_struc_annotations[position] = [position, 'none', solv_acc]

    
    # Finished looping through the file
    handle.close()
    
    # Save a table with the data
    df_sec_struc = pd.DataFrame.from_dict(sec_struc_annotations, orient='index', columns = ['Position', 'Secondary_structure', 'Solvent_accessibility'])
    df_sec_struc.to_csv(path_or_buf = os.path.join(out_folder, pdb_id + '_dssp_table.txt'), sep = '\t', index = False)
    
    ## Prepare the dictionary with the summary for each structure ##
    summary_dict[pdb_id] = {'Structure' : pdb_id, 
        'none' : 0, 'Missing' : 0, 'Alpha helix' : 0,
        'Beta bridge' : 0, 'Beta ladder' : 0, '3/10 helix' : 0, '5-helix' : 0, 'H-bonded turn' : 0,
        'Bend' : 0
    }
        
    ## Load the PDB file with the indicated regions
    pdb_interface_file = os.path.join(path_interfaces, pdb_id, 'dist_regions_' + pdb_id + '_bio_check_Repair.pdb')
    
    if os.path.exists(pdb_interface_file):
    
        dict_regions = {}

        ## Make a dictionary for the region of each position
        ## Loop through the lines
        with open(pdb_interface_file, 'r') as file_handle:
            # Parse pdb lines
            for line in file_handle:
                split_line = parse_pdb_line(line)
                # atom, atom_num, atom_name, resname, chain, res_num, x, y, z, bfactor

                res_num = int(split_line[5]) 
                region = float(split_line[9])

                # Save residue position and region
                dict_regions[res_num] = region
    
    summary_dict_interfaces[pdb_id] = {'Structure' : pdb_id, 
        'none' : 0, 'Missing' : 0, 'Alpha helix' : 0,
        'Beta bridge' : 0, 'Beta ladder' : 0, '3/10 helix' : 0, '5-helix' : 0, 'H-bonded turn' : 0,
        'Bend' : 0
    }
    
    # Loop through the residues and add to the count of each type of secondary structure
    for key, value_list in sec_struc_annotations.items():
        position = int(key)
        sec_struc_new = value_list[1]
        summary_dict[pdb_id][sec_struc_new] = summary_dict[pdb_id][sec_struc_new] + 1
    
        ## Use a second dictionary only for interface residues
        if os.path.exists(pdb_interface_file):

            if dict_regions.get(position, 0) in [0.75, 1.00]:
                summary_dict_interfaces[pdb_id][sec_struc_new] = summary_dict_interfaces[pdb_id][sec_struc_new] + 1
    


# Convert the summary dictionary to a data frame
# Save a table with the data
df_summary = pd.DataFrame.from_dict(summary_dict, orient='index')

df_summary

df_summary.to_csv(path_or_buf = os.path.join(out_folder, 'summary_table.txt'), sep = '\t', index = False)

# Convert the summary dictionary to a data frame
# Save a table with the data
df_summary_interfaces = pd.DataFrame.from_dict(summary_dict_interfaces, orient='index')
df_summary_interfaces.to_csv(path_or_buf = os.path.join(out_folder, 'summary_table_interfaces.txt'), sep = '\t', index = False)

