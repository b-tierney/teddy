#!/bin/bash

source activate /home/sez10/miniconda3_2/envs/r_env_glmnet

line=${1}

train_abundance_data=$(echo $line | awk -F ',' '{print $2}')
metadata_train=$(echo $line | awk -F ',' '{print $3}')

echo ${train_abundance_data}
echo ${metadata_train}
Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_ttests_2.R ${train_abundance_data} ${metadata_train}
