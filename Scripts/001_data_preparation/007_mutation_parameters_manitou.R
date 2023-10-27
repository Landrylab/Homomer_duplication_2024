########################################################
####            007_mutation_matrix_Manitou         ####
#### This script uses the gathered data to generate ####
#### the final matrices with mutational effects.    ####
########################################################

## This script must be run as follows:
## Rscript --vanilla <script_name> args1 args2 args3 args4

## This script will receive the following arguments:
# args1 = folder with the interface files (004_interfaces)
# args2 = folder with the gathered mutation files (006_gathered_mutations)
# args3 = pdb code to work with
# args4 = folder for the final matrices

args = commandArgs(trailingOnly = T)
interface_folder_path = args[1]
gathered_mut_path = args[2]
complex_name = args[3]
out_folder = args[4]

# Load libraries
library(ggplot2)
library(tidyverse)
library(magrittr)

# Create the output folder
if(!(dir.exists(out_folder))){
  dir.create(out_folder)
}
# Prepare a table of amino acid names
aa_three2one <- data.frame(cbind(c('A', 'R', 'D', 'N', 'C', 
                                   'E', 'Q', 'G', 'H', 'I', 
                                   'L', 'K', 'M', 'F', 'P',
                                   'S', 'T', 'W', 'Y', 'V'),
                                 c('ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                   'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                                   'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                   'SER', 'THR', 'TRP', 'TYR', 'VAL')))
colnames(aa_three2one) <- c('One-letter', 'Three-letter')

#### Define functions ####

# Write a function to format the data ####
preprocess_data <- function(HET_stab, HM_stab, HET_int, HM_int, complex_name){
  
  # Get only the needed columns
  HET_stab_small <- HET_stab %>% select(Chain, WT_res, Mut_res, Mean_ddG_stab_HET)
  HET_int_small <- HET_int %>% select(Chain, WT_res, Mut_res, Mean_ddG_int_HET)
  
  HM_stab_small <- HM_stab %>% select(Chain, WT_res, Mut_res, Mean_ddG_stab_HM)
  HM_int_small <- HM_int %>% select(Chain, WT_res, Mut_res, Mean_ddG_int_HM)
  
  # Join the two matrices
  tmp_mat <- inner_join(x = HET_stab_small, y = HM_stab_small, 
                        by = c("Chain" = "Chain", "WT_res" = "WT_res", "Mut_res" = "Mut_res"))
  tmp_mat2 <- inner_join(x = tmp_mat, y = HET_int_small,
                         by = c("Chain" = "Chain", "WT_res" = "WT_res", "Mut_res" = "Mut_res"))
  final_mat <- inner_join(x = tmp_mat2, y = HM_int_small,
                          by = c("Chain" = "Chain", "WT_res" = "WT_res", "Mut_res" = "Mut_res"))
  
  return(final_mat %>% distinct(WT_res, Mut_res, .keep_all = T))
}

# Write a function to prepare interface tables
get_interface_tables <- function(interface_data){
  interface_column <- strsplit(interface_data$Residues[1], split = ',|;')
  
  interface_table_tmp <- data.frame(interface_column[[1]])
  colnames(interface_table_tmp) <- c('Residue')
  
  interface_table_tmp %<>% rowwise() %>%
    mutate(Position = str_extract(Residue, pattern = '[0-9]+'),
           Three.letter = str_extract(Residue, pattern = '[A-Z]+'))
  
  interface_table <- left_join(interface_table_tmp, aa_three2one, by = c("Three.letter" = "Three-letter")) %>%
    unite(`One-letter`, Position, col = "WT_res", sep = '')
  
  return(interface_table)
  
}

#### Load data ####

## Interface files
interface_file <- paste(interface_folder_path, '/', complex_name,'/dist_regions_', complex_name, '_bio_check_Repair.pdb.txt', sep = '')
interface_data <- read_delim(interface_file, delim = '\t', col_names = c('Interface', 'Residues'))

## Stability ddGs
HET_stab<-read_delim(paste(gathered_mut_path, '/', complex_name, '_HET/stability_ddGs_tabs.txt', sep = ''),
                          delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_stab_HET', 'Std_dev_ddG_stab_HET', 'Min_ddG_stab_HET', 'Max_ddG_stab_HET'))
HM_stab<-read_delim(paste(gathered_mut_path, '/', complex_name, '_HM/stability_ddGs_tabs.txt', sep = ''),
                         delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_stab_HM', 'Std_dev_ddG_stab_HM', 'Min_ddG_stab_HM', 'Max_ddG_stab_HM'))

## Binding energy ddGs
HET_int<-read_delim(paste(gathered_mut_path, '/', complex_name, '_HET/A-B/interface_ddGs_tabs.txt', sep = ''), 
                         delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_int_HET', 'Std_dev_ddG_int_HET', 'Min_ddG_int_HET', 'Max_ddG_int_HET'))
HM_int<-read_delim(paste(gathered_mut_path, '/', complex_name, '_HM/A-B/interface_ddGs_tabs.txt', sep = ''), 
                        delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_int_HM', 'Std_dev_ddG_int_HM', 'Min_ddG_int_HM', 'Max_ddG_int_HM'))

# Build the final matrix
final_mat <- preprocess_data(HET_stab, HM_stab, HET_int, HM_int, complex_name)

# Save the data
write.table(final_mat, file = paste(out_folder, '/final_mat_', complex_name, '.txt', sep = ''), 
            quote = F, sep = '\t', row.names = F, col.names = T)






