#!/bin/bash

my_file=${1}

while read line
do
    abundance_data=$(echo ${line} | awk '{print $1}')
    train_metadata=$(echo ${line} | awk '{print $2}')
    test_metadata=$(echo ${line} | awk '{print $3}')
    mapping_data=$(echo ${line} | awk '{print $4}')
    training_subjects=$(echo ${line} | awk '{print $5}')
    testing_subjects=$(echo ${line} | awk '{print $6}')
    metadata_suffix=$(echo ${line} | awk '{print $7}')
    echo ${abundance_data}
    echo ${train_metadata}
    echo ${test_metadata}
    echo ${mapping_data}
    echo ${training_subjects}
    echo ${testing_subjects}
    /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_alignments/parsed_data_4/prep_abundance_for_voe_v4_bulk_training_data.bash ${abundance_data} ${train_metadata} ${test_metadata} ${mapping_data} ${training_subjects} ${testing_subjects} ${metadata_suffix}
done < ${my_file}
