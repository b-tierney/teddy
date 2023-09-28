#!/bin/bash

geneList=${1}
fastafile=${2}
output=${3}

source activate /home/sez10/miniconda3/envs/meta_assemblers

seqtk subseq ${fastafile} ${geneList} > ${output}
