#!/bin/bash

source activate r_env

abundance_data=${1}
metadata=${2}
mapping=${3}
metadata_suffix=${4}
Rscript prep_abundance_for_voe_v2.R ${abundance_data} ${metadata} ${mapping} ${metadata_suffix}
