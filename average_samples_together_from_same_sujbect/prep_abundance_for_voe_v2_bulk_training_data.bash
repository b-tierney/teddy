#!/bin/bash

source activate r_env

abundance_data=${1}
metadata=${2}
mapping=${3}
training_subjects=${4}
testing_sujbects=${5}
metadata_suffix=${6}
Rscript prep_abundance_for_voe_v2_bulk_training_data.R ${abundance_data} ${metadata} ${mapping} ${training_subjects} ${testing_sujbects} ${metadata_suffix}
