# # Symmetry check
# 
# A script to look at the distance between a residue in one chain of the homodimer and itself in the other subunit.

# Load libraries
import numpy as np
import pandas as pd
import matplotlib
import csv
import os
import sys
import subprocess
from Bio.PDB import PDBParser, PDBIO
import glob
from Bio import SeqIO
import re
from collections import OrderedDict
import math

## Some helper functions

# Use a class for PDB coordinates
class PDB_coordinates:
    def __init__(self, x_coord, y_coord, z_coord, residue, atomtype, bfactor, position, chain):
        '''The constructor for the PDB coordinates class'''
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
        self.residue = residue
        self.atomtype = atomtype
        self.bfactor = bfactor
        
        self.position = position
        self.chain = chain

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
                    # residue = atom_line[4] + atom_line[5] + atom_line[3]
                    
                    residue = atom_line[3]
                    position = atom_line[5]
                    chain = atom_line[4]
                    
                    
                    coord_dict[atom_line[4]].append(PDB_coordinates(x_coord, y_coord, z_coord, residue, 'CA', bfactor, position, chain))
            else:
                if coord_dict.get(atom_line[4], -1) == -1:
                        coord_dict[atom_line[4]] = []
                
                # Add the atom to the list
                x_coord = float(atom_line[6])
                y_coord = float(atom_line[7])
                z_coord = float(atom_line[8])
                bfactor = float(atom_line[9])
                atomtype = atom_line[2]
                # residue = atom_line[4] + atom_line[5] + atom_line[3]

                    
                residue = atom_line[3]
                position = atom_line[5]
                chain = atom_line[4]
                
                coord_dict[atom_line[4]].append(PDB_coordinates(x_coord, y_coord, z_coord, residue, atomtype, bfactor, position, chain))

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

## Folder for the PDB structures
path_structures = '../../Data/Structures/004_interfaces/'

list_structures = glob.glob(os.path.join(path_structures, '*'))

min_dist_list = []

## Loop through the list of structures
for in_folder in list_structures:
    ## Check the PDB ID and load the structure
    pdb_id = os.path.basename(in_folder)
    
    infile = os.path.join(in_folder, 'dist_regions_' + pdb_id + '_bio_check_Repair.pdb')
    
    if not os.path.exists(infile):
        continue
        
    ## Parse residues bfactor and positions
    struc_dict = parse_coordinates(infile, False)
    
    ## For residues at the interface, measure distance to their counterpart
    interface_dict = {'A':{}, 'B':{}}
    positions_interface = []
    for subunit, residue_data in struc_dict.items():
        for res in residue_data:
            ## Check if the residue is at the interface
            region = res.bfactor
            position = res.position
            residue = res.residue
            chain = res.chain
            
            if region >= 0.75:
                ## Save this residue to the list
                if not position in positions_interface:
                    positions_interface.append(position)
    
    ## Loop again to make sure every interface position is in the dictionary
    for subunit, residue_data in struc_dict.items():
        for res in residue_data:
            region = res.bfactor
            position = res.position
            residue = res.residue
            chain = res.chain
            
            if position in positions_interface:
                if interface_dict[chain].get(position, -1) == -1:
                    interface_dict[chain][position] = [res]
                else:
                    interface_dict[chain][position].append(res)
      
    ## Loop through the interface positions
    for position in positions_interface:
        
        # Initialize the minimum distance
        minimum_distance = 1000
        
        for atom1 in interface_dict['A'][position]:
            
            for atom2 in interface_dict['B'][position]:
                
                region = max(atom1.bfactor, atom2.bfactor)
                
                # Skip hydrogen atoms
                if not atom1.atomtype.startswith('H') and not atom2.atomtype.startswith('H'):
                    ## Measure the distance from atom1 to atom2
                    atom_distance = round(atom1.measure_distance(atom2), 2)
                
                    if atom_distance < minimum_distance:
                        minimum_distance = atom_distance
            
            
            ## End loop for atoms in chain B
        ## End loop for atoms in chain A
        
        ## Save the minimum distance between that residue and itself in the other subunit
        min_dist_list.append([pdb_id, atom1.position, atom1.residue, region, minimum_distance])
        
## Convert to a dataframe
min_dist_df = pd.DataFrame(min_dist_list, columns = ['PDB', 'Position', 'Residue', 'Region', 'Min_dist'])

min_dist_df.to_csv('../../Data/Structures/interface_distance_self.tsv', index = False, sep = '\t')

df_used_structures = pd.read_csv('../../Data/final_104_structures.txt', sep = '\t', names = ['PDB_ID'])
used_structures = [entry for entry in df_used_structures['PDB_ID']]
used_structures

list_check = [os.path.basename(entry) for entry in list_structures]

for entry in list_check:
    if not os.path.basename(entry) in used_structures:
        print(entry)

