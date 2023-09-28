#!/bin/bash

myfile=${1}
features=${2}
feature_to_keep=${3}
microbiome_features=${4}
thread_number=${5}
loss_function=${6}

while read line
do
/n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_lasso_binary_all_clinical.bash ${line} ${features} ${feature_to_keep} ${microbiome_features} ${thread_number} ${loss_function}
done < ${myfile}
