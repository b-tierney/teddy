#!/bin/bash

input_file=${1}
outputdir=${2}
diamond_index=${3}
diamond_fasta=${4}

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${input_file})

/n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/align_no_normalization.sh ${file} ${outputdir} ${diamond_index} ${diamond_fasta}
