#!/bin/bash

my_file=${1}

while read line
do
    abundance_data=$(echo ${line} | awk '{print $1}')
    metadata=$(echo ${line} | awk '{print $2}')
    mapping_data=$(echo ${line} | awk '{print $3}')
    training_subjects=$(echo ${line} | awk '{print $4}')
    testing_subjects=$(echo ${line} | awk '{print $5}')
    metadata_suffix=$(echo ${line} | awk '{print $6}')
    echo ${abundance_data}
    echo ${metadata}
    echo ${mapping_data}
    echo ${training_subjects}
    echo ${testing_subjects}
    ./prep_abundance_for_voe_v2_bulk_training_data.bash ${abundance_data} ${metadata} ${mapping_data} ${training_subjects} ${testing_subjects} ${metadata_suffix}
done < ${my_file}
