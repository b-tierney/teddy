#!/bin/bash

source activate r_env

mkdir -p ${2}
cd "${1}"

proportion_min_samples_abundant_in=${3}

#generate data locs and commands file
if [ ! -f datalocs ]
then
  ls healthy*_train_* | grep -v metadata > datalocs
fi
#while read p; do echo Rscript run_regressions.R $p healthy_*_1_metadata_filtered.rds; done<datalocs > commands

# get minimum value. will be saved as an rds file named "min_val"
#if [ ! -f min_val ]
#then
#  Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/get_minimum_value.R "${1}"
#fi

while read p
do
    Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_regressions_noVoE_rankNorm.R $p healthy*_train_*1_metadata_filtered.rds ${2} ${proportion_min_samples_abundant_in}
done < datalocs
