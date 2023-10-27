#!/bin/bash

########################################################
####	             008_simulation_workflow            ####
#### This script handles the simulations            ####
########################################################

#### Input parameters
bool_parametric=$1
main_result_folder=$2
path_final_matrices=$3
pdb_id=$4

mkdir -p ${main_result_folder}


#### Simulations with changes in synthesis rates ####
# Define the list of probabilities of mutations affecting expression we will use
# expression_probs_list=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
expression_probs_list=(0)

#### Simulations with different specific activities ####
## The two lists (aAAi_list and aABi_list) should have the same length

## Define the lists of values for specific activity
# aAAi_list=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1 1 1 1 1 1 1 1) ## Values used in figure 5E, 5F
aAAi_list=(1) ## Default value when specific activity is the same for homodimers and heterodimers

# aABi_list=(1 1 1 1 1 1 1 1 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2) ## Values used in figure 5E, 5F
aABi_list=(1) ## Default value when specific activity is the same for homodimers and heterodimers

#### Simulations with different starting values of folding energy and binding affinity ####
# The two lists (start_aff_list and start_stab_list) should have the same length

## Define the starting binding energy and affinity we will use
# start_aff_list=(-10 -10 -10 -5 -20) ## Values used for figure S6
# start_stab_list=(-5 -10 -2.5 -5 -5) ## Values used for figure S6

## Typical values for folding energy and binding affinity
start_aff_list=(-10)
start_stab_list=(-5)

#### Values for parametric simulations ####

## Define the values for the distribution of mutational effects for parametric distributions (means, shapes, correlations, sd)
## Mean values for the pooled distributions of mutational effects
mean_HET=0.2
mean_stab=2.6

## Values for the pooled distributions of mutational effects
sd_HET=1.2
sd_stab=4.6

## Shape for a multivariate normal distribution, 0 means no skew
shapes=(0 0 0)

# Correlations should be given in this order: (HMint_HETint HMint_stab HETint_stab)
correlations=(0.9 0.3 0.3)

mean_HM_list=(0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5) ## Values used in figures 2E, 4E, 5E

# std_dev_HM_list=(0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3) ## Values tested in Figure 2E
std_dev_HM_list=(2.4) ## Value for the pooled mutational effects of all structures

#### 008: Run simulations

## Loop through the provided values for specific activity
for i in $(seq 0 1 $((${#aAAi_list[@]} - 1)) )
do
	aABi="${aABi_list[$i]}"
	aAAi="${aAAi_list[$i]}"
	
	## Loop through the values for probability of changes in synthesis rates
	for j in $(seq 0 1 $((${#expression_probs_list[@]} - 1)) )
	do
        	expression_probs="${expression_probs_list[$j]}"
        	
        	## Loop through the values for starting binding affinity and starting folding energy
        	for k in $(seq 0 1 $((${#start_aff_list[@]} - 1)) )
        	do
            	start_aff="${start_aff_list[$k]}"
            	start_stab="${start_stab_list[$k]}"
	
            	## Check if we are running parametric simulations
            	if [ $bool_parametric -eq 1 ]
            	then
                	
                	## Loop through the values for the mean effect of mutations on HM binding
                	for l in $(seq 0 1 $((${#mean_HM_list[@]} - 1)))
                	do
                    	mean_HM="${mean_HM_list[$l]}"
                    	
                    	## Loop through the values for the standard deviations of mutations on HM binding
                    	for m in $(seq 0 1 $((${#std_dev_HM_list[@]} - 1)))
                    	do
                           sd_HM="${std_dev_HM_list[$m]}"
            	
                        	## Initialize the needed vectors for the parametric simulations
                        	means=($mean_HM $mean_HET $mean_stab)
                        	sd=(${sd_HM} ${sd_HET} ${sd_stab})
                        	
                        	out_dir_basename=aff${start_aff}_stab${start_stab}_intHM${mean_HM}_intHET${mean_HET}_stab${mean_stab}_exp${expression_probs}_sdHET${sd_HET}_sdHM${sd_HM}_HMact${aAAi}_HETact${aABi}

                           ## If this simulation has already been completed do not rerun it
                           if [ -f ${main_result_folder}/${out_dir_basename}/50/results.txt ]
                           then
                               echo 'Simulation ${out_dir_basename} already done, skipping'
                               echo -----
                               continue
                
                           fi

                           echo Submitting simulation ${out_dir_basename}

# Prepare the file with the input parameters
cat > input_params_simulations_slurm.py << end_file
## List of input parameters: 

## Indicate if we are doing parametric or single structure simulations
bool_parametric = $bool_parametric

# Alpha parameter (activity value for which the fitness function is maximized)
alpha = 60

# Beta parameter (fitness value when activity is twice or half of alpha)
beta = 0.5

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
aAAi = ${aAAi}

# Specific activity for heterodimers
aABi = ${aABi}

# Initial degradation rate for homodimers
dAAi = 1.3

# Probability of duplication
pdup = 1

# Initial value for binding energy
start_binding_energy = $start_aff

# Initial value for subunit stability
start_stab = $start_stab

# Number of replicates
# number_runs = 50
number_runs = 5

# Number of fixed mutations required to stop the simulation
# total_mutations = 200
total_mutations = 50

# Effective population size to use for the efficiency of selection
Ne = 500

# Output folder for the results
out_folder = '${main_result_folder}/${out_dir_basename}'

#### Parameters for gene expression changes #### 

# Probability of a mutation affecting gene expression
p_exp = $expression_probs

## Parameters for the distribution from which to sample the effects on gene expression
## based on data from Metzger et al., (2016). MBE

# Mean
mean_distr_exp = 0
# Standard deviation
sd_distr_exp = 0.025
# Skew
skew_distr_exp = -0.125

#### Parameters for distribution of mutational effects ####

## Need a covariance matrix between the three distributions (can calculate it from correlation coefficient + variances)
## Since r = cov(xy) / [sd(x) * sd *(y)]
## Then: cov(xy) = r * sd(x) * sd(y)

## Keep in mind I could get parameters for a skew normal distribution
## Individual distributions would become skew normal, covariance matrix would follow the same idea
## Just need to check how to sample from the multivariate skew normal distribution

## ddG HM int
mean_HM_int = ${means[0]}
sd_HM_int = ${sd[0]}
shape_HM_int = ${shapes[0]}

## ddG HET int
mean_HET_int = ${means[1]}
sd_HET_int = ${sd[1]}
shape_HET_int = ${shapes[1]}

## ddG stability
mean_stab = ${means[2]}
sd_stab = ${sd[2]}
shape_stab = ${shapes[2]}

## Correlations
r_HM_HET_int = ${correlations[0]}
r_HM_stab = ${correlations[1]}
r_HET_stab = ${correlations[2]}

end_file
        	
                          # Submit the simulations
                          ./008_simulations_slurm_manager.sh ${out_dir_basename} ${main_result_folder}/${out_dir_basename}
                          echo -----
    
                      done # End loop on standard deviation of effects of mutations on ddGbind,HM
                      
                  done # End loop on mean effect of mutations on ddGbind,HM
              else
              ## Work with single structure simulations now
              mkdir -p ${main_result_folder}/${pdb_id}
              out_dir_basename=${pdb_id}_aff${start_aff}_stab${start_stab}_exp${expression_probs}_HMact${aAAi}_HETact${aABi}

              ## If this simulation has already been completed do not rerun it
              if [ -f ${main_result_folder}/${pdb_id}/${out_dir_basename}/50/results.txt ]
              then
                  echo 'Simulation ${out_dir_basename} already done, skipping'
                  echo -----
                  continue
    
              fi
    
              echo Submitting simulation ${out_dir_basename}
    
# Prepare the file with the input parameters
cat > input_params_simulations_slurm.py << end_file
## List of input parameters: 

## Indicate if we are doing parametric or single structure simulations
bool_parametric = $bool_parametric

# Alpha parameter (activity value for which the fitness function is maximized)
alpha = 60

# Beta parameter (fitness value when activity is twice or half of alpha)
beta = 0.5

# Path to the matrix used as input
mut_matrix_file = '${path_final_matrices}/final_mat_${pdb_id}.txt'

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
aAAi = ${aAAi}

# Specific activity for heterodimers
aABi = ${aABi}

# Initial degradation rate for homodimers
dAAi = 1.3

# Probability of duplication
pdup = 1

# Initial value for binding energy
start_binding_energy = $start_aff

# Initial value for subunit stability
start_stab = $start_stab

# Number of replicates
# number_runs = 50
number_runs = 5

# Number of fixed mutations required to stop the simulation
# total_mutations = 200
total_mutations = 50

# Effective population size to use for the efficiency of selection
Ne = 500

# Output folder for the results
out_folder = '${main_result_folder}/${pdb_id}/${out_dir_basename}'

#### Parameters for gene expression changes #### 

# Probability of a mutation affecting gene expression
p_exp = $expression_probs

## Parameters for the distribution from which to sample the effects on gene expression
## based on data from Metzger et al., (2016). MBE

# Mean
mean_distr_exp = 0
# Standard deviation
sd_distr_exp = 0.025
# Skew
skew_distr_exp = -0.125


end_file
                
                # Submit the simulations
                ./008_simulations_slurm_manager.sh ${out_dir_basename} ${main_result_folder}/${pdb_id}/${out_dir_basename}
                echo -----
                        	
            	fi
        	
        	done # End loop on starting values for dGbind and dGstab
        	
    done # End loop on probabilities of mutations affecting gene expression
    
done # End loop on values for specific activity
	
	
	

	
