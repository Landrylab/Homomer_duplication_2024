# # Equation system solution space
# 
# This script will look at our system of equations and change the parameters to try to visualize how changes to each of them affect the solution space, that is, the proportion of complexes that are HET.
# 

# Import libraries
import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt
import seaborn as sns
import os
from datetime import date
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D

# Define helper functions
#### Prepare the functions that will apply the effects of the sampled mutations
def calculate_folded_fraction2(curr_stab):
    '''This function receives a starting deltaG of subunit stability. 
    With this, it calculates the new deltaG of stability and the fraction of
    folded proteins.
    '''
    folded_fraction = 1 / (1 + math.exp(curr_stab))

    return(folded_fraction)

def calculate_eq_constant2(curr_binding_energy, R, T):
    '''This function receives a starting binding energy.
    With this, it calculates the new equilibrium constant of association.
    '''    
    # Calculate the new equilibrium constant
    K = math.exp(-1 * curr_binding_energy/ (R * T))
    
    return(K)

# Define functions
def mean(lst):
    return(round(sum(lst) / len(lst), 2))

def dist_to_pos(number):
    if number.real >= 0:
        return abs(number.imag)
    else:
        return abs(number)

## Define some variables and constants
# Constants
R = 1.987e-3 # kcal / (mol * K)
T = 298 # K
dAA = 1.3
dBB = 1.3
dAB = 1.3
dA = 1.3
dB = 1.3

# Variables
sA = 100
sB = 100

curr_binding_energy_AA = -10
curr_binding_energy_BB = -10
curr_binding_energy_AB = -10

## Parameters for the fitness curve
alpha = 60
beta = 0.5

# Parameters for stability values
curr_stab_A=-5
curr_stab_B=-5

# Define the main function we will test
def equation_system(sA, sB, curr_stab_A, curr_stab_B, curr_binding_energy_AA, curr_binding_energy_BB, curr_binding_energy_AB):
    # Calculate the folded fractions
    folded_fraction_A = calculate_folded_fraction2(curr_stab_A)
    folded_fraction_B = calculate_folded_fraction2(curr_stab_B)
    
    # Calculate association constant
    kAA = calculate_eq_constant2(curr_binding_energy_AA, R, T)
    kBB = calculate_eq_constant2(curr_binding_energy_BB, R, T)
    kAB = calculate_eq_constant2(curr_binding_energy_AB, R, T)
    
    # Use the folded fraction to update synthesis rates
    sA = sA*(folded_fraction_A)
    sB = sB*(folded_fraction_B)

    # Double kAB to account for the increased collision probability
    kAB = kAB*2

    # Coefficients
    coeff4 = 2*dAA*kAA*(4*dAA*kAA*dBB*kBB/(dAB*kAB) - dAB*kAB)

    coeff3 = dA*(8*dAA*kAA*dBB*kBB/(dAB*kAB) - dAB*kAB) - 2*dB*dAA*kAA

    coeff2 = 2*dBB*kBB*(dA**2)/(dAB*kAB) + sA*(dAB*kAB - 8*dAA*kAA*dBB*kBB/(dAB*kAB)) - dA*dB - dAB*kAB*sB

    coeff1 = sA*(dB - 4*dA*dBB*kBB/(dAB*kAB))

    coeff0 = 2*(sA**2)*dBB*kBB/(dAB*kAB)

    coeffs = [coeff4, coeff3, coeff2, coeff1, coeff0]
    
    # Look for the roots
    cA = np.roots(coeffs)

    # There is thus a vector of 5 possible values for each concentration.

    cB = (sA/cA - 2*dAA*kAA*cA - dA)/(dAB*kAB)

    cAA = (cA**2)*kAA

    cBB = (cB**2)*kBB

    cAB = cA*cB*kAB
    
    # State which as the smallest sum of distances of concentrations
    # from the positive real numbers.
    most_positive_state = min(
        zip(cA, cAA, cAB, cBB, cB),
        key=lambda state: sum(map(dist_to_pos, state))
    )

    if not np.isclose(sum(map(dist_to_pos, most_positive_state)), 0):
        print("The numerical solving of the equilibrium state did not find positive concentrations.")
        return(0, 0, 0, 0, 0)

    # Take the absolute value of concentrations
    cA, cAA, cAB, cBB, cB = tuple(map(abs, most_positive_state))
    
    # Return the concentrations we need
    return(cA, cB, cAA, cBB, cAB)

equation_system(sA, sB, curr_stab_A, curr_stab_B, curr_binding_energy_AA, curr_binding_energy_BB, curr_binding_energy_AB)


# ## Test the effects of mutations on stability

stab_A_values = np.arange(-10, 10, 0.01)
stab_B_values = np.arange(-10, 10, 0.01)


list_results = []
for new_stab_A in stab_A_values:

    for new_stab_B in stab_B_values:

        # Feed the values into the system of equations
        cA, cB, cAA, cBB, cAB = equation_system(sA, sB, new_stab_A, new_stab_B, curr_binding_energy_AA, curr_binding_energy_BB, curr_binding_energy_AB)

        # Create a new line
        if (cA + cB + cAA + cAB + cBB) == 0:
            pct_AB = 0
        else:
            pct_AB = 100*cAB / (cA + cB + cAA + cAB + cBB)
            
        new_line = [new_stab_A, new_stab_B, cA, cB, cAA, cBB, cAB, pct_AB]
        
        list_results.append(new_line)

    
# Convert the list of lists to a dataframe
df = pd.DataFrame(list_results, columns = ['new_stab_A', 'new_stab_B', 'cA', 'cB', 'cAA', 'cBB', 'cAB', 'pct_AB'])

df.to_csv('../../Data/Results_solution_space/solutionSpace_subunitStab.tsv', 
         sep = '\t', index = False)


# ## Try with changes to the binding affinity of BB and AB

binding_AB_values = np.arange(-15 -5, 0.01)
binding_BB_values = np.arange(-15, -5, 0.01)

list_results = []
for new_binding_AB in binding_AB_values:

    for new_binding_BB in binding_BB_values:

        # Feed the values into the system of equations
        cA, cB, cAA, cBB, cAB = equation_system(sA, sB, curr_stab_A, curr_stab_B, curr_binding_energy_AA, new_binding_BB, new_binding_AB)

        # Create a new line
        if (cA + cB + cAA + cAB + cBB) == 0:
            pct_AB = 0
        else:
            pct_AB = 100*cAB / (cA + cB + cAA + cAB + cBB)
            
        new_line = [new_binding_AB, new_binding_BB, cA, cB, cAA, cBB, cAB, pct_AB]
        
        list_results.append(new_line)

    
# Convert the list of lists to a dataframe
df = pd.DataFrame(list_results, columns = ['new_binding_AB', 'new_binding_BB', 'cA', 'cB', 'cAA', 'cBB', 'cAB', 'pct_AB'])

df.to_csv('../../Data/Results_solution_space/solutionSpace_bindingEnergy.tsv', 
         sep = '\t', index = False)


# ## Panels for figures 5 and 6

# ## Modify synthesis rates on the x-axis (log2 scale) and binding affinities on the y-axis

## Synthesis rate of A will be constant
sA = 100
log2_values = list(np.arange(-2, 2, 0.05))

sB_values = [sA * (2**float(entry)) for entry in log2_values]

## Binding affinity of AA and BB will be constant
curr_binding_energy_AA = -10
curr_binding_energy_BB = -10
binding_values_AB = np.arange(-13, -7, 0.1)

list_results = []
# Loop through synthesis rate values
# for new_sB in sB_values:
for i in range(len(sB_values)):
    new_sB = sB_values[i]
    log2_val = log2_values[i]

    syn_ratio = round(new_sB / sA, 3)
    
    for new_binding_AB in binding_values_AB:

        binding_diff = new_binding_AB - curr_binding_energy_AA
        
        # Feed the values into the system of equations
        cA, cB, cAA, cBB, cAB = equation_system(sA, new_sB, curr_stab_A, curr_stab_B, curr_binding_energy_AA, curr_binding_energy_BB, new_binding_AB)

        # Create a new line
        pct_AB = 100*cAB / (cA + cB + cAA + cAB + cBB)
        new_line = [sA, new_sB, curr_binding_energy_AA, curr_binding_energy_BB, new_binding_AB,
                    cA, cB, cAA, cBB, cAB, pct_AB, syn_ratio, binding_diff, log2_val]
        
        list_results.append(new_line)

    
# Convert the list of lists to a dataframe
df = pd.DataFrame(list_results, columns = ['sA', 'sB', 'binding_energy_AA', 'binding_energy_BB', 'binding_energy_AB',
                                           'cA', 'cB', 'cAA', 'cBB', 'cAB', 'pct_AB', 'syn_ratio', 'binding_diff', 'log2_val'])

df

## Save the dataframe
df.to_csv('../../Data/Results_solution_space/solutionSpace_syn_ratio_diff_binding_log2.tsv', 
         sep = '\t', index = False)


# ## Modify specific activities on the x-axis and binding affinities on the y-axis

## Synthesis rate of A will be constant
sA = 100
sB = 100

hetBias_values = np.arange(-80, 81, 1)

## Binding affinity of AA and BB will be constant
curr_binding_energy_AA = -10
curr_binding_energy_BB = -10
binding_values_AB = np.arange(-13, -7, 0.1)

list_results = []
# Loop through synthesis rate values
for hetBias in hetBias_values:

    # Use the hetBias value to calculate specific activities of HM and HET
    low_act_value = 1 - abs(hetBias / 100)
    if hetBias < 0:
        # The HET should have lower activity
        aAA = 1
        aAB = low_act_value
    elif hetBias > 0:
        # The HET should have higher activity
        aAA = low_act_value
        aAB = 1
    elif hetBias == 0:
        # Both have the same activity
        aAA = 1
        aAB = 1
    
    for new_binding_AB in binding_values_AB:

        binding_diff = new_binding_AB - curr_binding_energy_AA
        
        # Feed the values into the system of equations
        cA, cB, cAA, cBB, cAB = equation_system(sA, sB, curr_stab_A, curr_stab_B, curr_binding_energy_AA, curr_binding_energy_BB, new_binding_AB)

        # Create a new line
        pct_AB = 100*cAB / (cA + cB + cAA + cAB + cBB)
        total_activity = cAA*aAA + cBB*aAA + cAB*aAB
        
        new_line = [sA, sB, curr_binding_energy_AA, curr_binding_energy_BB, new_binding_AB,
                    cA, cB, cAA, cBB, cAB, pct_AB, binding_diff, hetBias, total_activity]
        
        list_results.append(new_line)
    
# Convert the list of lists to a dataframe
df = pd.DataFrame(list_results, columns = ['sA', 'sB', 'binding_energy_AA', 'binding_energy_BB', 'binding_energy_AB',
                                           'cA', 'cB', 'cAA', 'cBB', 'cAB', 'pct_AB', 'binding_diff', 
                                          'hetBias', 'total_activity'])


df

## Save the dataframe
df.to_csv('../../Data/Results_solution_space/solutionSpace_hetBias_diff_binding.tsv', 
         sep = '\t', index = False)

