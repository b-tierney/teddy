#!/bin/bash

input_file=${1}
chr=${2}
ref_folder=${3}
output_prefix=${4}
module load java/jdk-1.8u112

# preprocess input data
myref=${ref_folder}/chr${chr}.1kg.phase3.v5a.vcf.gz
java -jar  /home/sez10/kostic_lab/software/conform-gt.24May16.cee.jar ref=${myref} gt=${input_file} chrom=${chr} out=${output_prefix}/${input_file%_vcf.vcf}_chr${chr}
