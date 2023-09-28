#!/bin/bash

abundance_data=${1}
metadata_file=${2}
HLA_list=${3}
condition_list=${4}
output_dir=${5}
abundance_type=${6}
source activate r_env

Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/humann_analysis/run_maaslin_species_pway.R ${abundance_data} ${metadata_file} ${HLA_list} ${condition_list} ${output_dir} ${abundance_type}
