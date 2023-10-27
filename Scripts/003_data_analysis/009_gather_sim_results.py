
# coding: utf-8

# # Gather simulation results
# 
# This script will look at the results from the simulations and gather them in a single file so that they are usable directly in R.
# 

# Load libraries
import numpy as np
import pandas as pd
import csv
import os
import glob
import re

## Define a path to the main folder for simulations that used PDB structures
main_folder = <path_to_folder>

### Define the name of the output file
out_file = os.path.join(main_folder, 'all_results_all_sims.tsv')
num_replicates = 50

bool_header = True

## Open the output file
with open(out_file, 'w') as out_file_handle:
    writer_handle = csv.writer(out_file_handle, delimiter = '\t')
    
    ## Loop through all the structures in the folder
    structure_list = glob.glob(os.path.join(main_folder, '*'))
    
    for structure_folder in structure_list:
        ## Save the structure ID
        structure_ID = os.path.basename(structure_folder)
    
        ## Loop through all the sets of conditions for that structure
        condition_list = glob.glob(os.path.join(structure_folder, '*'))
        for condition in condition_list:
            ## Save the current conditions
            condition_aff = re.search(pattern = 'aff(-[0-9]+)', string = os.path.basename(condition)).group(1)
            condition_stab = re.search(pattern = 'stab(-[0-9]+)', string = os.path.basename(condition)).group(1)
            
            ## Check if these are simulations that can modify expression values
            if re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)):
                exp_probs = re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)).group(1)
        
            ## Loop through all replicates
            for rep_number in range(1, 51):

                ## Load the data
                infile_path = os.path.join(condition, str(rep_number), 'results.txt')
                file_handle = csv.reader(open(infile_path, 'r'), delimiter = '\t')
                
                # Write the header only once
                if bool_header:
                    line = next(file_handle)
                    line.append('Replicate')
                    line.append('Complex')
                    line.append('Start_stab')
                    line.append('Start_aff')
                    
                    ## Add the column for probabilities of mutations affecting gene expression
                    if re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)):
                        line.append('Expression_mut_probs')
                    
                    writer_handle.writerow(line)
                    bool_header = False
                
                ## Loop through the lines in the file
                for line in file_handle:
                    ## Skip the header
                    if line[0] == 't':
                        continue
                    
                    ## Add columns for starting conditions and for the complex
                    line.append(rep_number)
                    line.append(structure_ID)
                    line.append(condition_stab)
                    line.append(condition_aff)
                    
                    ## Add the column for probabilities of mutations affecting gene expression
                    if re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)):
                        line.append(exp_probs)
                    
                    ## Write to file
                    writer_handle.writerow(line)
                
                    ## End loop for lines
                
                ## End loop for replicates
            
            ## End loop for conditions
            
        ## End loop for structures


# ## Use a separate loop to work with the simulations that used parametric distributions of mutational effects

main_folder = <path_to_main_folder>

condition_list = glob.glob(os.path.join(main_folder, 'aff*'))
condition_list

### Define the name of the output file
## Need to solve issue with permissions
out_file = os.path.join(main_folder, 'all_results_all_sims.tsv')
num_replicates = 50

bool_header = True

## Open the output file
with open(out_file, 'w') as out_file_handle:
    writer_handle = csv.writer(out_file_handle, delimiter = '\t')
   
    ## Loop through all the sets of conditions for that structure
    condition_list = glob.glob(os.path.join(main_folder, 'aff*'))
    for condition in condition_list:
        ## Save the current conditions
        condition_aff = re.search(pattern = 'aff(-[0-9]+)', string = os.path.basename(condition)).group(1)
        condition_stab = re.search(pattern = 'stab(-[0-9]+)', string = os.path.basename(condition)).group(1)

        ## Read the parameters used for the simulations
        intHM_param = re.search(pattern = 'intHM([0]?.[0-9]+)', string = os.path.basename(condition)).group(1)
        intHET_param = re.search(pattern = 'intHET([0-9]+.[0-9]+)', string = os.path.basename(condition)).group(1)
        stab_param = re.search(pattern = 'stab([0-9]+.[0-9]+)', string = os.path.basename(condition)).group(1)

        ## Check if these are simulations that can modify expression values
        if re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)):
            exp_probs = re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)).group(1)

        ## Add columns for standard deviation of mutations on the HM and HET
        if re.search(pattern = 'sdHET([0-9]?.[0-9]+)', string = os.path.basename(condition)):
            sdHET = re.search(pattern = 'sdHET([0-9]?.[0-9]+)', string = os.path.basename(condition)).group(1)

        if re.search(pattern = 'sdHM([0-9]?.[0-9]+)', string = os.path.basename(condition)):
            sdHM = re.search(pattern = 'sdHM([0-9]?.[0-9]+)', string = os.path.basename(condition)).group(1)
            
        ## Loop through all replicates
        list_subfolder = glob.glob(os.path.join(condition, '*'))
        list_replicates = []
        for entry in list_subfolder:
            basename_check = os.path.basename(entry)
            match_check = re.search(string = basename_check, pattern = '\d+')
            if match_check:
                if match_check.group(0) == basename_check:
                    list_replicates.append(int(basename_check))
        
        for rep_number in list_replicates:
            
            ## Load the data
            infile_path = os.path.join(condition, str(rep_number), 'results.txt')
            
            ## Skip folders for which the results file does not exist
            if not os.path.exists(infile_path):
                continue
            
            file_handle = csv.reader(open(infile_path, 'r'), delimiter = '\t')

            # Write the header only once
            if bool_header:
                line = next(file_handle)
                line.append('Replicate')
                # line.append('Complex')
                line.append('Start_stab')
                line.append('Start_aff')
                line.append('intHM_param')
                line.append('intHET_param')
                line.append('stab_param')
                
                ## Add the column for probabilities of mutations affecting gene expression
                if re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)):
                    line.append('Expression_mut_probs')
                
                ## Add columns for standard deviation of mutations on the HM and HET
                if re.search(pattern = 'sdHET([0-9]?.[0-9]+)', string = os.path.basename(condition)):
                    line.append('sdHET')
                    
                if re.search(pattern = 'sdHM([0-9]?.[0-9]+)', string = os.path.basename(condition)):
                    line.append('sdHM')
                
                writer_handle.writerow(line)
                bool_header = False

            ## Loop through the lines in the file
            for line in file_handle:
                ## Skip the header
                if line[0] == 't':
                    continue

                ## Add columns for starting conditions and for the complex
                line.append(rep_number)
                # line.append(structure_ID)
                line.append(condition_stab)
                line.append(condition_aff)
                line.append(intHM_param)
                line.append(intHET_param)
                line.append(stab_param)
                ## Add the column for probabilities of mutations affecting gene expression
                if re.search(pattern = 'exp([0-9].[0-9])', string = os.path.basename(condition)):
                    line.append(exp_probs)
                    
                ## Add columns for standard deviation of mutations on the HM and HET
                if re.search(pattern = 'sdHET([0-9]?.[0-9]+)', string = os.path.basename(condition)):
                    line.append(sdHET) ## Need to get these variables
                    
                if re.search(pattern = 'sdHM([0-9]?.[0-9]+)', string = os.path.basename(condition)):
                    line.append(sdHM)

                ## Write to file
                writer_handle.writerow(line)

                ## End loop for lines

            ## End loop for replicates

        ## End loop for conditions
            

