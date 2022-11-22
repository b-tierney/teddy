#!/bin/bash

input_prefix=${1}
ref_folder=${2}
map_input_folder=${3}
chr=${4}
start=${5}
stop=${6}
final_out_dir=${7}
thread_number=${8}
module load java/jdk-1.8u112

myref=${ref_folder}/chr${chr}.1kg.phase3.v5a.bref3
# now phase and impute
java -Xmx50g -jar /home/sez10/kostic_lab/software/beagle.28Jun21.220.jar gt=${input_prefix}${chr}.vcf.gz ref=${myref} out=${final_out_dir}/chr${chr}_${start}_${stop}_imputed map=${map_input_folder}/plink.chr${chr}.GRCh37.map chrom=${chr}:${start}-${stop} impute=true nthreads=${thread_number}
