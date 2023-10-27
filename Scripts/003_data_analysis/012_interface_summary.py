# # Interface summary
# 
# This script will look at the data from the files that called interfaces to make tables with:
# - Numbers of core residues
# - Numbers of rim residues
# - Percentage of residues that are in the core
# - Percentage of residues that are in the rim
# 

# Load libraries
import numpy as np
import pandas as pd
import csv
import os
import glob
import re

import math
from collections import OrderedDict
from Bio.PDB import *

## Load some helper functions

# Use a class for PDB coordinates
class PDB_coordinates:
    def __init__(self, x_coord, y_coord, z_coord, residue, atomtype, bfactor):
        '''The constructor for the PDB coordinates class'''
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
        self.residue = residue
        self.atomtype = atomtype
        self.bfactor = bfactor

    def measure_distance(self, target):
        '''A function that measures the distance between two PDB_coordinates objects.'''
        x_dist = (self.x - target.x)**2
        y_dist = (self.y - target.y)**2
        z_dist = (self.z - target.z)**2
        return math.sqrt(x_dist + y_dist + z_dist)

##################################

def parse_coordinates(infile, alpha_only):
    '''This function will parse a pdb file using the PDB_coordinates class.
    If the "alpha_only" argument is set to True the function only reads alpha carbons ("CA").
    If it is set to False it will read all atoms.
    '''
    handle = open(infile, 'r')
    coord_dict = OrderedDict()

    for line in handle:
        if line.startswith('ATOM'):

            atom_line = parse_pdb_line(line)

            if alpha_only:

                if atom_line[2] == 'CA':

                    if coord_dict.get(atom_line[4], -1) == -1:
                        coord_dict[atom_line[4]] = []
                    
                    # Add the atom to the list
                    x_coord = float(atom_line[6])
                    y_coord = float(atom_line[7])
                    z_coord = float(atom_line[8])
                    bfactor = float(atom_line[9])
                    residue = atom_line[4] + atom_line[5] + atom_line[3]
                    coord_dict[atom_line[4]].append(PDB_coordinates(x_coord, y_coord, z_coord, residue, 'CA', bfactor))
            else:
                if coord_dict.get(atom_line[4], -1) == -1:
                        coord_dict[atom_line[4]] = []
                
                # Add the atom to the list
                x_coord = float(atom_line[6])
                y_coord = float(atom_line[7])
                z_coord = float(atom_line[8])
                bfactor = float(atom_line[9])
                atomtype = atom_line[2]
                residue = atom_line[4] + atom_line[5] + atom_line[3]
                coord_dict[atom_line[4]].append(PDB_coordinates(x_coord, y_coord, z_coord, residue, atomtype, bfactor))

    return coord_dict

##################################


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
    bfactor = pdb_line[60:66].strip(' ')

    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z, bfactor]

# Go to the main folder
main_folder = '../../Data/Structures/004_interfaces'

list_files = glob.glob(os.path.join(main_folder, '*'))

list_files

out_list = []

# Loop through each of the structures
for pdb_folder in list_files:

    pdb_id = os.path.basename(pdb_folder)

    # Parse the file
    pdb_file = os.path.join(pdb_folder, 'dist_regions_' + pdb_id + '_bio_check_Repair.pdb')
    if os.path.isfile(pdb_file):
        pdb_coords = parse_coordinates(pdb_file, True)

        # Count number of residues per subunit (should be the same)
        length_struc = len(pdb_coords['A'])

        core_res = 0
        rim_res = 0
        # Count number of residues at the interface or the rim
        for residue in pdb_coords['A']:
            if residue.bfactor == 0.75:
                rim_res += 1
            elif residue.bfactor == 1:
                core_res += 1

        out_list.append([pdb_id, length_struc, core_res, rim_res])

df = pd.DataFrame(out_list, columns = ['Structure', 'Total_residues', 'Interface_residues', 'Rim_residues'])

df

df.to_csv('../../Data/Structures/interface_summary.tsv', 
         sep = '\t',  header = True, index = False)

# ## Write a dataframe with the stickiness values for residues at the interface

## Load the stickiness scale
aa_properties = pd.read_csv('../../Data/Levy2012_propensity.tsv', sep = '\t')
aa_properties

codontable_standard = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

## Will also need the reverse dictionary
aa_three2one = {
    'ALA': 'A', 'ILE': 'I', 'MET': 'M', 'THR': 'T', 'ASN': 'N',
    'LYS': 'K', 'SER': 'S', 'ARG': 'R', 'LEU': 'L', 'PRO': 'P',
    'HIS': 'H', 'GLN': 'Q', 'VAL': 'V', 'ASP': 'D', 'GLU': 'E',
    'GLY': 'G', 'PHE': 'F', 'TYR': 'Y', 'CYS': 'C', 'TRP': 'W'
}

out_list = []

# Loop through each of the structures
for pdb_folder in list_files:

    pdb_id = os.path.basename(pdb_folder)

    # Parse the file
    pdb_file = os.path.join(pdb_folder, 'dist_regions_' + pdb_id + '_bio_check_Repair.pdb')
    if os.path.isfile(pdb_file):
        pdb_coords = parse_coordinates(pdb_file, True)

        # Count number of residues at the interface or the rim
        for residue_i in pdb_coords['A']:
            # Extract the residue type and the position
            res_matches = re.search(pattern = '(\d+)(.*)', string = residue_i.residue)
            position = res_matches.group(1)
            res_type = aa_three2one[res_matches.group(2)]
            
            ## Check the stickiness value for that residue
            aa_values = aa_properties[aa_properties['Aminoacid.1.letter'] == res_type]
            for index, line in aa_values.iterrows():
                stickiness_value = line['levy_propensity']
            
            if residue_i.bfactor == 0.75:
                region = 'Rim'
                out_list.append([pdb_id, position, res_type, stickiness_value, region])
                
                
            elif residue_i.bfactor == 1:
                region = 'Core'
                out_list.append([pdb_id, position, res_type, stickiness_value, region])

out_df = pd.DataFrame(out_list, columns = ['PDB', 'Position', 'Res_type', 'Levy_propensity', 'Region'])
out_df

means_stickiness = out_df.groupby('PDB')['Levy_propensity'].mean()
means_stickiness.describe()

out_df.to_csv('../../Data/Structures/interface_stickiness.tsv', 
         sep = '\t',  header = True, index = False)

