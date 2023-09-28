#!/bin/bash

input_file=${1}
proportion_filter=${2}
thread_number=${3}

source activate /home/sez10/miniconda3_2/envs/r_env_glmnet

while read line
do
  echo ${line}
  Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/filter_abundance_metadata_data_v2_species_pathways.R ${line} ${proportion_filter} ${thread_number}
done < ${input_file}

#conda deactivate

#input_folder_basename=$(basename ${input_foler})
#train_input=$(echo ${input_foler}"/"${input_folder_basename}"filtered_abundance_train.csv")
#test_input=$(echo ${input_foler}"/"${input_folder_basename}"filtered_abundance_test.csv")

#train_input_names=$(echo ${input_foler}"/"${input_folder_basename}"filtered_abundance_train_names.txt")
#test_input_names=$(echo ${input_foler}"/"${input_folder_basename}"filtered_abundance_test_names.txt")


#train_output=$(echo ${input_foler}"/"${input_folder_basename}"filtered_transformed_abundance_train.csv")
#test_output=$(echo ${input_foler}"/"${input_folder_basename}"filtered_transformed_abundance_test.csv")
#source activate python_env 

#echo $train_input
#echo $test_input
#echo $train_output
#echo $test_output

#python /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/normalize_abundance_data_v2.py ${train_input} ${test_input} ${train_input_names} ${test_input_names} ${train_output} ${test_output}
