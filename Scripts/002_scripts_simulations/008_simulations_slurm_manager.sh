#!/bin/bash

## This code will manage the script for simulations by preparing a job and submitting it to SLURM
## The only parameters it will receive are a tag for the job name and the path to the output folder

job_tag=$1
out_dir=$2

# Define the SLURM queue
slurm_queue=<slurm_queue>

mkdir -p $out_dir
cp 008_simulations_slurm.py $out_dir
cp input_params_simulations_slurm.py $out_dir

cd $out_dir

cat > ${job_tag}.sbatch << EOF
#!/bin/bash

#SBATCH -J ${job_tag}
#SBATCH -o ${job_tag}.out
#SBATCH -c 1
#SBATCH -p $slurm_queue
#SBATCH --time=1-00:00
#SBATCH --mem=2000

# Load Python and libraries
module load python/3.8.10
module load scipy-stack/2020a

# Run the simulations
python3.8 008_simulations_slurm.py

EOF

sbatch ${job_tag}.sbatch

