#!/bin/bash

my_file=${1}

while read line
do
    abundance_data=$(echo ${line} | awk '{print $1}')
    metadata=$(echo ${line} | awk '{print $2}')
    mapping_data=$(echo ${line} | awk '{print $3}')
    metadata_suffix=$(echo ${line} | awk '{print $4}')
    echo ${abundance_data}
    echo ${metadata}
    echo ${mapping_data}
    ./prep_abundance_for_voe_v2.bash ${abundance_data} ${metadata} ${mapping_data} ${metadata_suffix}
done < ${my_file}
