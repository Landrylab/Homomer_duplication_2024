#!/bin/bash

# This code receives the name of a pdb file
# and passes it on to the FoldX RepairPDB function.
# The SLURM queue parameter should be modified as needed
# $1 = path to the input PDB without the 'pdb' extension
# $2 = path to the output directory
# $3 = number of times to repeat 

full_path=$1
repeat_num=$3

prot=$(basename $1)

mkdir -p $2

cat > ${prot}_repair.sbatch << EOF
#!/bin/bash

#SBATCH -D $PWD
#SBATCH -J ${prot}_repair
#SBATCH -o ${prot}_repair.out
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --time=1-00:00
#SBATCH --mem=51200

cp ${full_path}.pdb $2

cd $2

### Start the loop

for i in \$(seq 1 ${repeat_num})
do
    # Make sure that the rotabase.txt file needed by FoldX can be linked to this folder
    ln -s `which rotabase.txt` rotabase.txt
    
    FoldX --command=RepairPDB --pdb=${prot}.pdb --ionStrength=0.05 --pH=7 --water=CRYSTAL --vdwDesign=2 --out-pdb=true --pdbHydrogens=false > ${prot}_Repair.log
    
    mkdir ${prot}_Repair_\${i}
    
    mv ${prot}_Repair.* ${prot}_Repair_\${i}
    
    rm ${prot}.pdb
    
    cp ${prot}_Repair_\${i}/${prot}_Repair.pdb ./${prot}.pdb

done

EOF

sbatch ${prot}_repair.sbatch

