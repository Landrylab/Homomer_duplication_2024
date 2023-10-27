#!/bin/bash

#######################################################
####            006_gather_mutations               ####
#### This script will receive a folder of the      ####
#### output of mutagenesis runs to organize all    ####
#### the results for analysis.                     ####
#######################################################

#### Input arguments:
# $1 = absolute path to the folder with the mutatex results
# $2 = absolute path to the output folder

mkdir -p $2

cd $1

MAINDIR=$PWD

# Read the mutation list
mut_list_path=${MAINDIR}/mutation_list.txt

#### Start with stability effects ####

# Go into the folder that contains all the data
cd results/mutation_ddgs
res_folder=$(ls | grep -v 'final_averages')
cd $res_folder

# Concatenate all the files, removing the headers
# Add the columns for the starting residue (repeat file names, 20 times each)
# and the new residue (repeat mutation list)
for infile in $(ls)
do
    
    # Get the residue name
    old_resname=$(echo $infile | cut -f 1 -d '_' | cut -c 1,3-)
    chain=$(echo $infile | cut -f 1 -d '_' | cut -c 2)
    for i in $(seq 1 1 20)
    do
        echo -e $chain'\t'$old_resname >> column_old_resname.txt
    done
    
    # Get the data for these mutations without the headers
    tail -n +2 $infile > tmp_file.txt
    
    # Add new lines to the output file
    # Column names will be:
    # Old residue
    # New residue
    # Average ddG
    # Std. dev. of ddG (5 replicates)
    # Minimum ddG
    # Maximum ddG
    paste column_old_resname.txt ${mut_list_path} tmp_file.txt >> $2/stability_ddGs.txt
    
    rm column_old_resname.txt tmp_file.txt
    
done

sed 's/ /\t/g' $2/stability_ddGs.txt > $2/stability_ddGs_tabs.txt

#### Start with interfaces ####

# Go into the folder that contains all the data
cd ${MAINDIR}/results/interface_ddgs
res_folder=$(ls | grep -v 'final_averages')
cd $res_folder

# Loop through each of the interfaces
for interface_folder in $(ls)
do
    cd $interface_folder
    mkdir $2/${interface_folder}
    
    # Concatenate all the files, removing the headers
    # Add the columns for the starting residue (repeat file names, 20 times each)
    # and the new residue (repeat mutation list)

    for infile in $(ls)
    do
    
        # Get the residue name
        old_resname=$(echo $infile | cut -f 1 -d '_' | cut -c 1,3-)
        chain=$(echo $infile | cut -f 1 -d '_' | cut -c 2)
        for i in $(seq 1 1 20)
        do
            echo -e $chain'\t'$old_resname >> column_old_resname.txt
        done
        
        # Get the data for these mutations without the headers
        tail -n +2 $infile > tmp_file.txt
        
        # Add new lines to the output file
        # Column names will be:
        # Old residue
        # New residue
        # Average ddG
        # Std. dev. of ddG (5 replicates)
        # Minimum ddG
        # Maximum ddG
        paste column_old_resname.txt $mut_list_path tmp_file.txt >> $2/${interface_folder}/interface_ddGs.txt
        
        rm column_old_resname.txt tmp_file.txt
        
    done

    sed 's/ /\t/g' $2/${interface_folder}/interface_ddGs.txt > $2/${interface_folder}/interface_ddGs_tabs.txt

    cd ..
done