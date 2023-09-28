#!/bin/bash

metadata_train=${1}
metadata_test=${2}
feature_list=${3}

Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_logistic_regression.R ${metadata_train} ${metadata_test} ${feature_list}
