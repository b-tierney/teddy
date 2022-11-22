#!/bin/bash

input_file=${1}
out_dir=${2}
adapter_file=${3}
bitmask_file=${4}
srprism_file=${5}
tmp_dir=${6}
threads=${7}
dmd_index=${8}
fasta_ref=${9}

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ${input_file})

/n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/diabimmune_align_and_preprocess.bash ${file} ${out_dir} ${adapter_file} ${bitmask_file} ${srprism_file} ${tmp_dir} ${threads} ${dmd_index} ${fasta_ref}
