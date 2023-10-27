
#######################################################
####            008_simulations_code               ####
#### This script will receive a matrix of          ####
#### mutational effects, starting parameters, and  ####
#### an output directory to run simulations        ####
#######################################################

# Import libraries
import numpy as np
import pandas as pd
import math
import os
import shutil
from scipy.stats import skewnorm
from   scipy.stats import (multivariate_normal as mvn,
                           norm)
from   scipy.stats._multivariate import _squeeze_output

# Define the class for the multivariate skewnormal distributions
# Reference: http://gregorygundersen.com/blog/2020/12/29/multivariate-skew-normal/
class multivariate_skewnorm:
    
    def __init__(self, shape, cov=None):
        self.dim   = len(shape)
        self.shape = np.asarray(shape)
        self.mean  = np.zeros(self.dim)
        self.cov   = np.eye(self.dim) if cov is None else np.asarray(cov)

    def pdf(self, x):
        return np.exp(self.logpdf(x))
        
    def logpdf(self, x):
        x    = mvn._process_quantiles(x, self.dim)
        pdf  = mvn(self.mean, self.cov).logpdf(x)
        cdf  = norm(0, 1).logcdf(np.dot(x, self.shape))
        return _squeeze_output(np.log(2) + pdf + cdf)
    
    def rvs_fast(self, size=1):
        aCa      = self.shape @ self.cov @ self.shape
        delta    = (1 / np.sqrt(1 + aCa)) * self.cov @ self.shape
        cov_star = np.block([[np.ones(1),     delta],
                             [delta[:, None], self.cov]])
        x        = mvn(np.zeros(self.dim+1), cov_star).rvs(size)
        x0, x1   = x[:, 0], x[:, 1:]
        inds     = x0 <= 0
        x1[inds] = -1 * x1[inds]
        return x1
    
#### Run the script that defines the input parameters
from input_params_simulations_slurm import *

#### Define functions ####

# Define functions
def mean(lst):
    return(round(sum(lst) / len(lst), 2))
    
#### Prepare the functions that will apply the effects of the sampled mutations
def calculate_folded_fraction(curr_stab):
    '''This function receives a starting deltaG of subunit stability. 
    With this, it calculates the new deltaG of stability and the fraction of
    folded proteins.
    '''
    folded_fraction = 1 / (1 + math.exp(curr_stab))

    return(folded_fraction)

def calculate_eq_constant(curr_binding_energy, R, T):
    '''This function receives a starting binding energy.
    With this, it calculates the new equilibrium constant of association.
    '''    
    # Calculate the new equilibrium constant
    K = math.exp(-1 * curr_binding_energy/ (R * T))
    
    return(K)
    
# Function that computes the distance from a complex number
# to the closest positive real number.
def dist_to_pos(number):
    if number.real >= 0:
        return abs(number.imag)
    else:
        return abs(number)

# Function to compute monomer and dimer abundances, their specific
# activities and the log2 of fitness. Only works before duplication.
def predup_phenotype(sA, dA, aA, dAA, curr_binding_energy, curr_stab, aAA):

    # Calculate association constant
    kAA = calculate_eq_constant(curr_binding_energy, R, T)

    # Calculate the folded fraction
    folded_fraction = calculate_folded_fraction(curr_stab)
    # Use the folded fraction to update sA
    sA = sA*(folded_fraction)

    # Equilibrium concentration of monomer A
    cA = (-dA + (dA**2 + 8*dAA*kAA*sA)**0.5)/(4*dAA*kAA)

    # Equilibrium concentration of dimer AA
    cAA = (cA**2)*kAA

    # Total activity of monomers and dimers
    activity = cA*aA + cAA*aAA

    # log2 of fitness
    logw = np.log2(beta**(np.log2(activity/alpha)**2))

    return (cA, aA, cAA, aAA, logw, activity)

## Same as above, but only works after duplication.
# Params will now contain: (17 terms)
# sA, dA, aA, dAA, curr_binding_energy_AA, curr_stab_A, aAA,
# sB, dB, aB, dBB, curr_binding_energy_BB, curr_stab_B, aBB,
# dAB, curr_binding_energy_AB, aAB
def postdup_phenotype(sA, dA, aA, dAA, curr_binding_energy_AA, curr_stab_A, aAA, sB, dB, aB, dBB, curr_binding_energy_BB, curr_stab_B, aBB, dAB, curr_binding_energy_AB, aAB):      

    # Calculate the folded fractions
    folded_fraction_A = calculate_folded_fraction(curr_stab_A)
    folded_fraction_B = calculate_folded_fraction(curr_stab_B)

    # Use the folded fraction to update synthesis rates
    sA = sA*(folded_fraction_A)
    sB = sB*(folded_fraction_B)

    # Calculate association constant
    kAA = calculate_eq_constant(curr_binding_energy_AA, R, T)
    kBB = calculate_eq_constant(curr_binding_energy_BB, R, T)
    kAB = calculate_eq_constant(curr_binding_energy_AB, R, T)

    # Double kAB to account for the increased collision probability
    kAB = kAB*2

    # Coefficients of the polynomial to solve numerically

    coeff4 = 2*dAA*kAA*(4*dAA*kAA*dBB*kBB/(dAB*kAB) - dAB*kAB)

    coeff3 = dA*(8*dAA*kAA*dBB*kBB/(dAB*kAB) - dAB*kAB) - 2*dB*dAA*kAA

    coeff2 = 2*dBB*kBB*(dA**2)/(dAB*kAB) + sA*(dAB*kAB - 8*dAA*kAA*dBB*kBB/(dAB*kAB)) - dA*dB - dAB*kAB*sB

    coeff1 = sA*(dB - 4*dA*dBB*kBB/(dAB*kAB))

    coeff0 = 2*(sA**2)*dBB*kBB/(dAB*kAB)


    coeffs = [coeff4, coeff3, coeff2, coeff1, coeff0]


    # The roots of the polynomial are a vector of 5 possible values
    # of the concentration of the monomer A.

    cA = np.roots(coeffs)

    # There is thus a vector of 5 possible values for each concentration.

    cB = (sA/cA - 2*dAA*kAA*cA - dA)/(dAB*kAB)

    cAA = (cA**2)*kAA

    cBB = (cB**2)*kBB

    cAB = cA*cB*kAB

    # Among these 5 possible states, we have to figure out which one
    # is made of positive real numbers (within a small margin of error)

    # State which as the smallest sum of distances of concentrations
    # from the positive real numbers.
    most_positive_state = min(
        zip(cA, cAA, cAB, cBB, cB),
        key=lambda state: sum(map(dist_to_pos, state))
    )

    if not np.isclose(sum(map(dist_to_pos, most_positive_state)), 0):
        # raise Exception("The numerical solving of the equilibrium state did not find positive concentrations.")
        cA, cAA, cAB, cBB, cB = tuple(map(abs, most_positive_state))
        return(cA, cAA, cAB, cBB, cB, aA, aAA, aAB, aBB, aB, -100, 0)

    # Take the absolute value of concentrations
    cA, cAA, cAB, cBB, cB = tuple(map(abs, most_positive_state))

    # Total activity of monomers and dimers
    activity = cA*aA + cAA*aAA + cAB*aAB + cBB*aBB + cB*aB

    # log2 of fitness
    logw = np.log2(beta**(np.log2(activity/alpha)**2))

    return (cA, cAA, cAB, cBB, cB, aA, aAA, aAB, aBB, aB, logw, activity)


#### Start of evolve function ####
    
def evolve(sAi, dAi, aAi, aAAi, dAAi, aABi, pdup,
           Ne, start_binding_energy, start_stab,
           p_exp, mean_distr_exp, sd_distr_exp, skew_distr_exp,
           bool_parametric,
           *args, **kwargs):
    
    # sAi: initial synthesis rate of monomer A (units of concentration/time)
    # dAi: degradation rate of monomer A (units of time**-1)
    # aAi: specific activity of monomer A (units of activity/concentration)
    # dAAi: degradation rate of dimer AA (units of time**-1)
    # aAAi: specific activity of dimer AA (units of activity/concentration)
    # aABi: specific activity of dimer AB (units of activity/concentration)
    
    # pdup: probability of each mutation to be a duplication
    # Ne: effective population size
    # start_binding_energy: the starting value for binding energy
    # start_stab: the starting value for stability
    
    ## Arguments for the multivariate skew normal distribution of mutational effects
    # vector_means: the vector of means
    # vector_shapes: the shapes (skew)
    # cov_matrix: covariance matrix to indicate variances and correlation
    
    ## Arguments to sample effects on gene expression
    # p_exp: the probability of mutations affecting expression instead of the coding sequence
    # mean_distr_exp: the mean of the distribution from which to sample expression effects
    # sd_distr_exp: the standard deviation of the distribution from which to sample expression effects
    # skew_distr_exp: the skew of the distribution from which to sample expression effects
    
    ## Other arguments for parametric simulations:
    # vector_means: vector of means of the multivariate normal distribution of effects (ddGbind_HM, ddGbind_HET, ddGfold)
    # vector_shapes = vector of shapes of the multivariate normal distribution of effects (ddGbind_HM, ddGbind_HET, ddGfold)
    # cov_matrix = covariance matrix of the multivariate normal distribution of effects
    
    ## Other arguments for single structure simulations:
    # mutation_matrix = table with the distribution of mutational effects for a structure, from which mutations will be sampled
    
    ## Get arguments for parametric simulations if provided
    vector_means = kwargs.get('vector_means', None)
    vector_shapes = kwargs.get('vector_shapes', None)
    cov_matrix = kwargs.get('cov_matrix', None)
    
    ## Get arguments for single structure simulations if provided
    mutation_matrix = kwargs.get('mutation_matrix', None)
    
    #### Initialize the multivariate skew normal distribution to sample mutational effects for parametric simulations
    if bool_parametric:
        dist_mut_eff = multivariate_skewnorm(shape = vector_shapes, cov = cov_matrix)
    
    # A counter to keep track of how many mutations have fixed so far
    counter = 0
    
    # Initialize the variables for current binding energy and stability as the starting values
    curr_binding_energy = start_binding_energy
    curr_stab = start_stab
    
    ## Array of mutated parameters
    params = [sAi, dAi, aAi, dAAi, curr_binding_energy, curr_stab, aAAi]
    
    # Initial phenotype
    cA, aA, cAA, aAA, logw, activity = predup_phenotype(*params)
    
    # Initially, there is no most recently fixed mutation, so there is
    # no selection coefficient and no synonymous status of the mutation
    s = None
    synon = None
    
    # The initial time is 0
    t = 0
    
    # The gene has not fixed a duplication yet
    dup = False

    fixation = False

    # A bool to check if a mutation on exp fixed
    bool_exp = False
    
    # Each turn of this loop is the fixation of a mutation. It ends when
    # the fixed mutation is a duplication
    while not dup:
        
        # Yield the state of the simulation now that a new mutation
        # has fixed (or the initial state if t==0).
        yield dict(
            params=params, cA=cA, cAA=cAA, cAB=None, cBB=None,
            cB=None, aA=aA, aAA=aAA, aAB=None, aBB=None, aB=None,
            logw=logw, s=s, synon=synon, genecopy=None, t=t,
            activity=activity
        )
        
        # Check if the previous fixed mutation affected expression
        if fixation and bool_exp:
            print('Mutation on expression fixed')
        
        fixation = False
        
        # A bool to check if a mutation on exp fixed
        bool_exp = False
        
        # Each turn of this loop is a suggested mutation. The loop ends
        # when a mutation fixes.
        while not fixation:
            
            # Each suggested mutation counts as one unit of time,
            # regardless of whether it fixes or not
            t += 1
            
            # With probability pdup, the suggested mutation is a
            # duplication.
            dup = np.random.random() < pdup
            
            # Initialize new_params as the current params before applying changes
            new_params = params
            
            # This if-else block creates a new_log_params array
            # representing the phenotype of the suggested mutant.
            if dup:

                # If the mutation is a duplication, the rate of protein
                # synthesis doubles
                new_params = new_params * np.array([2, 1, 1, 1, 1, 1, 1])
                
                synon = None  
                
                # The duplication does not affect folding or affinity, so I will save the ddG = 0
                ddG_stab = 0
                ddG_binding_energy = 0
                
                ## Mutant phenotype
                # Prepare an array of parameter changes (sAi, aAi, dAAi, binding_energy, stability, aAAi)
                params_change = np.array([0, 0, 0, 0, ddG_binding_energy, ddG_stab, 0])

            else:
                
                #### Add an if to check if the mutation will affect gene expression or the other parameters
                if np.random.random() < p_exp:
                    
                    bool_exp = True
                    
                    # Sample from skew normal distribution
                    # Arguments are: mean (loc), skew (a), sd (scale), number of samples (size)
                    mult_mut_eff_exp = skewnorm.rvs(loc = mean_distr_exp, a = skew_distr_exp, scale = sd_distr_exp, size = 1)
                    
                    # Calculate the change in expression by multiplying the above effect times
                    # the previous value of expression
                    mut_eff_exp = mult_mut_eff_exp * new_params[0]
                    
                    # Prepare an array of parameter changes (sAi, aAi, dAAi, binding_energy, stability, aAAi)
                    params_change = np.array([mut_eff_exp, 0, 0, 0, 0, 0, 0])
                
                else:
                    ## Sample effects on the coding sequence
                
                    if bool_parametric:
                        ## Sample from the parametric distribution of mutational effects:
                        
                        # The sampled effects will be [HM_int, HET_int, stab].
                        sampled_row = dist_mut_eff.rvs_fast(size = 2)
                        
                        ddG_stab = float(sampled_row[0, 2] + vector_means[2])
                        ddG_binding_energy = float(sampled_row[0, 0] + vector_means[0])
                        
                    else:
                        ## Sample from the single structure distribution of mutational effects:
                        sampled_row = mutation_matrix.sample()
                
                        ddG_stab = float(sampled_row['Monomer_ddG'])
                        ddG_binding_energy = float(sampled_row['Mean_ddG_int_HM'])                        
                
                    ## Mutant phenotype
                    # Prepare an array of parameter changes (sAi, aAi, dAAi, binding_energy, stability, aAAi)
                    params_change = np.array([0, 0, 0, 0, ddG_binding_energy, ddG_stab, 0])
            
            # Save the new params (update dAi, binding_energy, stability or expression as needed)
            new_params = new_params + params_change
            
            # Replace this line to get the new phenotype
            # new_cA, new_aA, new_cAA, new_aAA, new_logw = predup_phenotype(*np.exp2(new_log_params))
            new_cA, new_aA, new_cAA, new_aAA, new_logw, new_activity = predup_phenotype(*new_params)
            
            # Selection coefficient
            s = np.exp2(new_logw - logw) - 1
            
            # Probability of fixation
            if s == 0:
                pfix = 1/Ne
            else:
                pfix = (1 - np.exp(-2*s))/(1 - np.exp(-2*s*Ne))
            
            # Whether or not the mutation fixes
            fixation = np.random.random() <= pfix
        
        # Now that a mutation has fixed, update log parameters
        params = new_params
        counter = counter + 1
        print('Mutation accepted', counter)
        
        # Update phenotype
        cA, aA, cAA, aAA, logw, activity = new_cA, new_aA, new_cAA, new_aAA, new_logw, new_activity
    
    ### A duplication has now fixed
    
    print('-------')
    print('Duplication fixed')
    print('-------')
    
    # Now that a duplication has fixed, we complexify log_params
    
    # Divide synthesis rate by 2 to recover the pre-duplication array
    params[0] = params[0] / 2
    
    ## Extend log_params to include parameters for one more monomer, one
    ## more homodimer and the heterodimer    
    ## Params will now contain: (17 terms)
    # sA, dA, aA, dAA, curr_binding_energy_AA, curr_stab_A, aAA,
    # sB, dB, aB, dBB, curr_binding_energy_BB, curr_stab_B, aBB,
    # dAB, curr_binding_energy_AB, aAB
    params = np.concatenate([params, params, params[[3, 4]], np.array([aABi])])
    
    # Update the phenotype
    cA, cAA, cAB, cBB, cB, aA, aAA, aAB, aBB, aB, logw, activity = postdup_phenotype(*params)
    
    genecopy = None
    
    # Each turn of this endless loop is the fixation of a point mutation
    # in one copy of the gene.
    while True:
        
        # Yield the state of the simulation now that a new mutation
        # has fixed.
        yield dict(
            params=params, cA=cA, cAA=cAA, cAB=cAB, cBB=cBB,
            cB=cB, aA=aA, aAA=aAA, aAB=aAB, aBB=aBB, aB=aB, logw=logw,
            s=s, synon=synon, genecopy=genecopy, t=t,
            activity=activity
        )
        
        # Check if the previous fixed mutation affected expression
        if fixation and bool_exp:
            print('Mutation on expression fixed')
        
        fixation = False
        
        # A bool to check if a mutation on exp fixed
        bool_exp = False
        
        # Each turn of this loop is a suggested mutation. The loop ends
        # when a mutation fixes.
        while not fixation:
            
            # Each suggested mutation counts as 0.5 unit of time,
            # regardless of whether it fixes or not (0.5 is to account
            # for the doubling of the mutation rate after duplication)
            t += 0.5
            
            # With probability 0.5, the suggested point mutation
            # is in copy A of the gene, otherwise in copy B.
            genecopy = np.random.choice(["A", "B"])
            
            #### Add an if to check if the mutation will affect gene expression or the other parameters
            if np.random.random() < p_exp:

                bool_exp = True
                
                # Sample from skew normal distribution
                # Arguments are: mean (loc), skew (a), sd (scale), number of samples (size)
                mult_mut_eff_exp = skewnorm.rvs(loc = mean_distr_exp, a = skew_distr_exp, scale = sd_distr_exp, size = 1)

                if genecopy == "A":

                    # Calculate the change in expression by multiplying the above effect times
                    # the previous value of expression
                    mut_eff_exp = mult_mut_eff_exp * params[0]
                    params_change = np.array([mut_eff_exp, 0, 0, 0, 0,0, 0,
                                             0, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0])
                else:

                    # Calculate the change in expression by multiplying the above effect times
                    # the previous value of expression
                    mut_eff_exp = mult_mut_eff_exp * params[7]
                    params_change = np.array([0, 0, 0, 0, 0,0, 0,
                                             mut_eff_exp, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0])
            else:
                
            
                bool_exp = False
                
                # Sample from the distribution of mutational effects
                # This will give me:
                # - ddG for stability
                # - ddG for binding energy
                
                if bool_parametric:
                    ## Sample from the parametric distribution of mutational effects:
                    # The sampled effects will be [HM_int, HET_int, stab].
                    sampled_row = dist_mut_eff.rvs_fast(size = 2)
                        
                    # Sampled mutation for A
                    ddG_binding_energy_AA = float(sampled_row[0, 0] + vector_means[0])
                    ddG_binding_energy_A_HET = float(sampled_row[0, 1] + vector_means[1])
                    ddG_stab_A = float(sampled_row[0, 2] + vector_means[2])
                    
                    # Sampled mutation for B
                    ddG_binding_energy_BB = float(sampled_row[1, 0] + vector_means[0])
                    ddG_binding_energy_B_HET = float(sampled_row[1, 1] + vector_means[1])
                    ddG_stab_B = float(sampled_row[1, 2] + vector_means[2])
                    
                else:
                    ## Sample from the single structure distribution of mutational effects:
                    sampled_row_A = mutation_matrix.sample() 
                    ddG_stab_A = float(sampled_row_A['Monomer_ddG'])
                    ddG_binding_energy_AA = float(sampled_row_A['Mean_ddG_int_HM'])
                    ddG_binding_energy_A_HET = float(sampled_row_A['Mean_ddG_int_HET'])
        
                    sampled_row_B = mutation_matrix.sample() 
                    ddG_stab_B = float(sampled_row_B['Monomer_ddG'])
                    ddG_binding_energy_BB = float(sampled_row_B['Mean_ddG_int_HM'])
                    ddG_binding_energy_B_HET = float(sampled_row_B['Mean_ddG_int_HET'])
                
                ## Defne an array of mutational effects
                ## Params contains: (17 terms)
                # sA, dA, aA, dAA, curr_binding_energy_AA, curr_stab_A, aAA,
                # sB, dB, aB, dBB, curr_binding_energy_BB, curr_stab_B, aBB,
                # dAB, curr_binding_energy_AB, aAB
                if genecopy == "A":
                    params_change = np.array([0, 0, 0, 0, ddG_binding_energy_AA, ddG_stab_A, 0,
                                             0, 0, 0, 0, 0, 0, 0,
                                             0, ddG_binding_energy_A_HET, 0])
                
                else:
                    params_change = np.array([0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, ddG_binding_energy_BB, ddG_stab_B, 0,
                             0, ddG_binding_energy_B_HET, 0])

            new_params = params + params_change
            
            # Mutant phenotype
            new_cA, new_cAA, new_cAB, new_cBB, new_cB, new_aA, new_aAA, new_aAB, new_aBB, new_aB, new_logw, new_activity = postdup_phenotype(*new_params)
            
            # Handle sets of parameters for which the system has no solution
            if new_logw:

                # Selection coefficient
                s = np.exp2(new_logw - logw) - 1

                # Probability of fixation
                if s == 0:
                    pfix = 1/Ne
                else:
                    pfix = (1 - np.exp(-2*s))/(1 - np.exp(-2*s*Ne))

                # Whether or not the mutation fixes
                fixation = np.random.random() <= pfix
                
            else:
                # If there was no solution, just discard the mutation and try with a new one
                pfix = 0
                fixation = np.random.random() <= pfix
        
        # Since the mutation as fixed, update log parameters
        params = new_params
        counter = counter + 1
        print('Mutation accepted', counter, genecopy)
        
        # Update phenotype
        cA, cAA, cAB, cBB, cB, aA, aAA, aAB, aBB, aB, logw, activity = new_cA, new_cAA, new_cAB, new_cBB, new_cB, new_aA, new_aAA, new_aAB, new_aBB, new_aB, new_logw, new_activity

#### End of evolve function ####

def write_results(results, results_folder):
    '''This function will organize the results into a pandas dataframe and write them to the output folder.'''
    os.makedirs(results_folder, exist_ok = True)
    
    # Get the concentrations
    t = [d["t"] for d in results]
    cAA = [d["cAA"] for d in results]
    cA = [d["cA"] for d in results]
    cBB = [d["cBB"] for d in results]
    cB = [d["cB"] for d in results]
    cAB = [d["cAB"] for d in results]
    
    # Get the rest of the parameters
    ## Params contains 7 terms:
    # sAi, dAi, aAi, dAAi, curr_binding_energy, curr_stab, aAAi
    
    ## Params will now contain: (17 terms)
    # sA, dA, aA, dAA, curr_binding_energy_AA, curr_stab_A, aAA,
    # sB, dB, aB, dBB, curr_binding_energy_BB, curr_stab_B, aBB,
    # dAB, curr_binding_energy_AB, aAB
    params_list = [d["params"] for d in results]
    new_params_list = []
    
    # Loop through the params list
    for i in range(len(params_list)):
        # Add the value of time
        new_params = np.insert(np.array(params_list[i]), 0, t[i])
    
        # If this sublist has 8 elements (including time), add 'None' ten times at the end to respect the number of columns
        if len(new_params) == 8:
            new_params = np.append(new_params, [None, None, None, None, None, None, None, None, None, None])
    
        new_params_list.append(new_params)
        
    # Concatenate the entries of the list as rows to build the dataframe
    df_results = pd.DataFrame(np.array(new_params_list))
    
    # Set column names
    df_results.columns = ['t', 'sA', 'dA', 'aA', 'dAA', 'curr_binding_energy_AA', 
                         'curr_stab_A', 'aAA','sB', 'dB', 'aB', 'dBB', 
                         'curr_binding_energy_BB', 'curr_stab_B', 'aBB','dAB', 'curr_binding_energy_AB', 'aAB']
    
    # Add concentration columns at the end
    df_results['cAA'] = cAA
    df_results['cA'] = cA
    df_results['cBB'] = cBB
    df_results['cB'] = cB
    df_results['cAB'] = cAB
    
    # Write dataframe
    results_file = os.path.join(results_folder, "results.txt")
    df_results.to_csv(results_file, sep = '\t', header = True, index = False)

#### Run simulations ####

if bool_parametric:    
    ## Prepare inputs for multivariate normal distributions (parametric distributions)
    
    vector_means = [mean_HM_int, mean_HET_int, mean_stab]
    vector_corr = [r_HM_HET_int, r_HM_stab, r_HET_stab]
    vector_shapes = [shape_HM_int, shape_HET_int, shape_stab]
    
    # Prepare covariance matrix
    cov_matrix = np.array([(0.0, 0.0, 0.0), 
                          (0.0, 0.0, 0.0),
                          (0.0, 0.0, 0.0)])
    
    cov_matrix[0,0] = sd_HM_int**2
    cov_matrix[0,1] = r_HM_HET_int * sd_HM_int * sd_HET_int
    cov_matrix[0,2] = r_HM_stab * sd_HM_int * sd_stab
    
    cov_matrix[1,0] = r_HM_HET_int * sd_HM_int * sd_HET_int
    cov_matrix[1,1] = sd_HET_int**2
    cov_matrix[1,2] = r_HET_stab * sd_HET_int * sd_stab
    
    cov_matrix[2,0] = r_HM_stab * sd_HM_int * sd_stab
    cov_matrix[2,1] = r_HET_stab * sd_HET_int * sd_stab
    cov_matrix[2,2] = sd_stab**2

else:
    #### Load the mutation matrix for the structure we are working with
    
    # Load matrix
    mutation_matrix = pd.read_csv(mut_matrix_file, delimiter = '\t')
    mutation_matrix['Monomer_ddG'] = mutation_matrix['Mean_ddG_stab_HET']


# Create output folder and save parameters there
os.makedirs(out_folder, exist_ok = True)

# Copy the file with the input parameters to the output folder
readme_handle = open(os.path.join(out_folder, 'README.txt'), 'w')

readme_handle.write('Set of starting parameters used for these runs:\n')
readme_handle.write('sAi = ' + str(sAi) + '\n')
readme_handle.write('dAi = ' + str(dAi) + '\n')
readme_handle.write('aAi = ' + str(aAi)+ '\n') 
readme_handle.write('aAAi = ' + str(aAAi)+ '\n')
readme_handle.write('dAAi = ' + str(dAAi)+ '\n')
readme_handle.write('pdup = ' + str(pdup)+ '\n')
readme_handle.write('Ne = ' + str(Ne)+ '\n')
readme_handle.write('start_binding_energy = ' + str(start_binding_energy)+ '\n')
readme_handle.write('start_stab = ' + str(start_stab)+ '\n')
readme_handle.write('number_runs = ' + str(number_runs)+ '\n')
readme_handle.write('fixed_mutations_per_run = ' + str(total_mutations)+ '\n')

if bool_parametric:
    readme_handle.write('vector_means = ' + str(vector_means)+ '\n')
    readme_handle.write('vector_var = ' + str(cov_matrix.diagonal())+ '\n')
    readme_handle.write('vector_shapes = ' + str(vector_shapes)+ '\n')
    readme_handle.write('vector_corr = ' + str(vector_corr)+ '\n')
else:
    readme_handle.write('mutation_matrix = ' + mut_matrix_file+ '\n')


readme_handle.close()

for i in range(1, number_runs + 1):
    
    ## Run parametric simulations
    if bool_parametric:
    
        gen = evolve(sAi = sAi,
                 dAi = dAi,
                 aAi = aAi, 
                 aAAi = aAAi,
                 dAAi= dAAi,
                 aABi = aABi,
                 pdup = pdup,
                 Ne= Ne,
                 start_binding_energy = start_binding_energy, 
                 start_stab = start_stab,
                 p_exp = p_exp,
                 mean_distr_exp = mean_distr_exp, 
                 sd_distr_exp = sd_distr_exp, 
                 skew_distr_exp = skew_distr_exp,
                 bool_parametric = bool_parametric,
                 
                 vector_means = vector_means,
                 vector_shapes = vector_shapes,
                 cov_matrix = cov_matrix,
                 )
    ## Run simulations with the mutational effects for a single structure
    else:
        gen = evolve(sAi = sAi,
                 dAi = dAi,
                 aAi = aAi, 
                 aAAi = aAAi,
                 dAAi= dAAi,
                 aABi = aABi,
                 pdup = pdup,
                 Ne= Ne,
                 start_binding_energy = start_binding_energy, 
                 start_stab = start_stab,
                 p_exp = p_exp,
                 mean_distr_exp = mean_distr_exp, 
                 sd_distr_exp = sd_distr_exp, 
                 skew_distr_exp = skew_distr_exp,
                 bool_parametric = bool_parametric,
                 mutation_matrix = mutation_matrix
                 )


    results = []
    # Rerun this line to extend the simulation by N fixed mutations
    results.extend(next(gen) for j in range(total_mutations))
    
        
    # Create output folder for this run and save results
    new_folder = os.path.join(out_folder, str(i))
    write_results(results, new_folder)


