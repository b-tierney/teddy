#!/bin/bash

source activate r_env_glmnet

line=${1}

test_abundance_data=$(echo $line | awk -F ',' '{print $1}')
train_abundance_data=$(echo $line | awk -F ',' '{print $2}')
metadata_train=$(echo $line | awk -F ',' '{print $3}')
metadata_test=$(echo $line | awk -F ',' '{print $4}')
train_sujbects1=$(echo $line | awk -F ',' '{print $5}')
train_sujbects2=$(echo $line | awk -F ',' '{print $6}')
train_sujbects3=$(echo $line | awk -F ',' '{print $7}')
test_subjects=$(echo $line | awk -F ',' '{print $8}')
feature_list=${2} # example would be microbiome,number_autoantibodies,fdr,grs2
features_to_keep=${3} # example number_autoantibodies,fdr,grs2
microbiome_features=${4} # either a file name or the word all
number_threads=${5}
loss_function=${6}

echo ${test_abundance_data}
echo ${train_abundance_data}

Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_rf_survival.R ${test_abundance_data} ${train_abundance_data} ${metadata_train} ${metadata_test} ${train_sujbects1} ${train_sujbects2} ${train_sujbects3} ${test_subjects} ${feature_list} ${features_to_keep} ${microbiome_features} ${number_threads} ${loss_function}
