#!/bin/bash

source activate /home/sez10/miniconda3_2/envs/r_env

abundance_data=${1}
train_metadata=${2}
test_metadata=${3}
mapping=${4}
training_subjects=${5}
testing_sujbects=${6}
metadata_suffix=${7}
Rscript /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/parsed_data_4/prep_abundance_for_voe_v4_bulk_training_data.R ${abundance_data} ${train_metadata} ${test_metadata} ${mapping} ${training_subjects} ${testing_sujbects} ${metadata_suffix}
