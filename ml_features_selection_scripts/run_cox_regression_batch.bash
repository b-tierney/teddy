#!/bin/bash                                                                                                                                                                                                                                 
source activate r_env_glmnet

input_file=${1}
feature_list=${2} # example would be number_autoantibodies,fdr,grs2                                                                                                                                                             

while read line
do
  metadata_train=$(echo $line | awk -F ',' '{print $3}')
  metadata_test=$(echo $line | awk -F ',' '{print $4}')
  echo $metadata_train
  echo $metadata_test
  Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_cox_regression.R ${metadata_train} ${metadata_test} ${feature_list}
done < ${input_file}
