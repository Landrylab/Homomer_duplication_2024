## List of input parameters: 

## Indicate if we are doing parametric or single structure simulations
bool_parametric = 0

# Alpha parameter (activity value for which the fitness function is maximized)
alpha = 60

# Beta parameter (fitness value when activity is twice or half of alpha)
beta = 0.5

# Path to the matrix used as input
mut_matrix_file = '<path_to_mutational_matrix>

# Gas constant
R = 1.987e-3

# Temperature
T = 298

# Initial synthesis rate
sAi = 100

# Initial degradation rate for monomers
dAi = 1.3

# Initial specific activity for monomers
aAi = 0.1

# Initial specific activity for homodimers
aAAi = 1

# Specific activity for heterodimers
aABi = 0.5

# Initial degradation rate for homodimers
dAAi = 1.3

# Probability of duplication
pdup = 1

# Initial value for binding energy
start_binding_energy = -10

# Initial value for subunit stability
start_stab = -5

# Number of replicates
# number_runs = 50
number_runs = 5

# Number of fixed mutations required to stop the simulation
# total_mutations = 200
total_mutations = 50

# Effective population size to use for the efficiency of selection
Ne = 500

# Output folder for the results
out_folder = <path_to_output_folder>

#### Parameters for gene expression changes #### 

# Probability of a mutation affecting gene expression
p_exp = 0

## Parameters for the distribution from which to sample the effects on gene expression
## based on data from Metzger et al., (2016). MBE

# Mean
mean_distr_exp = 0
# Standard deviation
sd_distr_exp = 0.025
# Skew
skew_distr_exp = -0.125


