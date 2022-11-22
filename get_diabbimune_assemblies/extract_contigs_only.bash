#!/bin/bash

input=${1}
output_dir=${2}

mkdir -p ${output_dir}
#download file
wget ${input}

# get the name of the file to decompress
sampleName=$(echo ${input##*/})
# decompress
fna_name=$(tar -xzvf ${sampleName}  --wildcards --no-anchored '*fna' -C ./)

# remove decompressed file
rm ${sampleName}

path_name=$(dirname ${fna_name})
parent_dir=$(basename ${path_name})
# remove suffix _prokka_out
sample_only=${parent_dir%_prokka_out}

mv ${fna_name} ${output_dir}/${sample_only}.fna
