#!/bin/bash

ped_file=${1}
input_dir=${2}
source activate eigensoft
cd ${input_dir}
smartpca -p ${ped_file}
