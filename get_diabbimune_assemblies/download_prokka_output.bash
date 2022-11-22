#!/bin/bash

input=${1}
output_dir=${2}

mkdir -p ${output_dir}
#download file
wget ${input}

# get the name of the file to decompress
sampleName=$(echo ${input##*/})
# decompress
tar_output=$(tar -xvzf ${sampleName})
# remove decompressed file
rm ${sampleName}
path_of_file=$(echo ${tar_output} | cut -f 1 -d ' ')

parent_dir=$(basename ${path_of_file})
# remove suffix _prokka_out
sample_only=${parent_dir%_prokka_out}

mkdir -p ${output_dir}/${parent_dir}

mv ${path_of_file}/* ${output_dir}/${parent_dir}
# prokka output has random file names. lets name them after sample instead
for f in ${output_dir}/${parent_dir}/*
do
# first remove path
myextension=$(echo ${f##*.})
mv $f ${output_dir}/${parent_dir}/${sample_only}.${myextension}
done
