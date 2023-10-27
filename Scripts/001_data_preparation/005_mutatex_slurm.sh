#!/bin/bash

# This code will receive an input PDB file for mutagenesis and the option
# to do it considering multimers or not. Input arguments are:
# $1 = path to the input PDB (without the pdb extension)
# $2 = path to the output directory
# $3 = multimer or not (0 for no multimers, 1 for multimers)
# $4 = number of cores to use

# Organize the output directory
mkdir -p $2
cp $1.pdb $2

cp mutation_list.txt $2
cd $2

## Define paths to mutatex home directory and the required virtual environment
MUTATEX_HOME=<path_to_mutatex_home>
MUTATEX_VENV=<path_to_mutatex_venv>

# Define SLURM queue
slurm_queue=<slurm_queue>

prot=$(basename $1)

# Set the flag for multimers or no multimers
if [ $3 -eq 0 ]
then
    mut_arg='--no-multimers'
fi

cat > mut_$3_${prot}.sbatch << EOF
#!/bin/bash

#SBATCH -D $PWD
#SBATCH -J mut_$3_${prot}
#SBATCH -o mut_$3_${prot}.out
#SBATCH -c $4
#SBATCH -p $slurm_queue
#SBATCH --time=1-00:00
#SBATCH --mem=8000

unset $PYTHONPATH
source $MUTATEX_VENV

mutatex ${prot}.pdb \
        -p $4 \
        -m mutation_list.txt \
        -x $FOLDX_BINARY \
        -f suite5 \
        -R ${MUTATEX_HOME}/templates/foldxsuite4/repair_runfile_template.txt \
        -M ${MUTATEX_HOME}/templates/foldxsuite4/mutate_runfile_template.txt \
        -I ${MUTATEX_HOME}/templates/foldxsuite4/interface_runfile_template.txt \
        -B $mut_arg


EOF

sbatch mut_$3_${prot}.sbatch

