# # Gather mutation matrices
# 
# This script will look at the matrices of mutational effects and gather them in a single file so that they are usable directly in R

# Load libraries
import numpy as np
import pandas as pd
import csv
import os
import glob
import re

## Define a path to the main folder
main_folder = <path_to_main_folder>

### Define the name of the output file
## Need to solve issue with permissions
out_file = os.path.join(main_folder, 'all_mut_matrices.tsv')

header = ['Chain', 'WT_res', 'Position', 'Mut_res', 'Mean_ddG_stab_HET', 
          'Mean_ddG_stab_HM', 'Mean_ddG_int_HET', 'Mean_ddG_int_HM', 'Complex']

## Open the output file
with open(out_file, 'w') as out_file_handle:
    writer_handle = csv.writer(out_file_handle, delimiter = '\t')
    
    writer_handle.writerow(header)
    
    ## Loop through all the mutational matrices
    structure_list = glob.glob(os.path.join(main_folder, 'final_mat_*'))
    
    for structure_file in structure_list:
        ## Save the structure ID
        structure_ID = re.search(pattern = '_([a-z0-9]+).txt', string = os.path.basename(structure_file)).group(1)
    
        ## Load the data
        # infile_path = os.path.join(condition, str(rep_number), 'results.txt')
        file_handle = csv.reader(open(structure_file, 'r'), delimiter = '\t')

        ## Loop through the lines in the file
        for line in file_handle:
            ## Skip the header
            if line[0] == 'Chain':
                continue

            ## Distinguish the three cases I have
            if structure_ID == '1p6o':
                ## Add columns for starting conditions and for the complex
                line.append(structure_ID)

                ## Write to file
                writer_handle.writerow(line)
            else:
                ## Separate the WT residue and the position
                WT_res = line[1][0]
                mut_pos = line[1][1:]
                
                line = [line[0], WT_res, mut_pos] + line[2:]
                
                if structure_ID == '1m38' or structure_ID == '4fgw':
                    # Remove the final column
                    line = line[0:-1]
                
                ## Add columns for starting conditions and for the complex
                line.append(structure_ID)

                ## Write to file
                writer_handle.writerow(line)


            ## End loop for lines
                    
        ## End loop for structures

