#!/bin/bash

vcf_file=${1}
bref_file=${2}
module load java/jdk-1.8u112
java -jar  /home/sez10/kostic_lab/software/bref3.28Jun21.220.jar ${vcf_file} > ${bref_file}
