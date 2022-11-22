#!/bin/bash
source activate /home/sez10/miniconda3/envs/meta_assemblers

dirname=${1}
output_folder=${2}
adapters_reference=${3}
bitmask_file=${4}
sprism_file=${5}
temp_dir=${6}
thread_number=${7}
mkdir -p ${output_folder}/${dirname}
#pairedEnd=0

export OMP_NUM_THREADS=${thread_number}

module load sratoolkit/2.10.7

prefetch --force no --ngc /n/scratch3/users/b/btt5/TEDDY/prj_20109.ngc -O ${output_folder}/${dirname} $dirname
fastq-dump --split-files --ngc /n/scratch3/users/b/btt5/TEDDY/prj_20109.ngc --outdir ${output_folder}/${dirname} ${output_folder}/${dirname}/${dirname}_dbGaP-20109.sra
# remove sra file
rm ${output_folder}/${dirname}/${dirname}_dbGaP-20109.sra

f1=${output_folder}/${dirname}/${dirname}_dbGaP-20109_1.fastq
f2=${output_folder}/${dirname}/${dirname}_dbGaP-20109_2.fastq

bbduk.sh -Xmx8g in1=${f1} in2=${f2} out1=${output_folder}/${dirname}/${dirname}_1.trimmed.fastq out2=${output_folder}/${dirname}/${dirname}_2.trimmed.fastq ref=${adapters_reference} threads=${thread_number} qin=33 ktrim=r k=23 mink=11 hdist=1 tpe tbo
# remove raw reads
rm ${f1}
rm ${f2}
# remove contamination
bmtagger.sh -b ${bitmask_file} -x ${sprism_file} -T ${temp_dir} -q1 -1 ${output_folder}/${dirname}/${dirname}_1.trimmed.fastq -2 ${output_folder}/${dirname}/${dirname}_2.trimmed.fastq -o ${output_folder}/${dirname}/${dirname}_human_free -X
rm ${output_folder}/${dirname}/${dirname}_1.trimmed.fastq
rm ${output_folder}/${dirname}/${dirname}_2.trimmed.fastq

## run megahit 
time megahit --mem-flag 2 --num-cpu-threads ${thread_number} -1 ${output_folder}/${dirname}/${dirname}_human_free_1.fastq -2 ${output_folder}/${dirname}/${dirname}_human_free_2.fastq -o ${output_folder}/${dirname}/${dirname}_mh_out
## cat the 2 files together for aligning to diamond later
cat ${output_folder}/${dirname}/${dirname}_human_free_1.fastq ${output_folder}/${dirname}/${dirname}_human_free_2.fastq > ${output_folder}/${dirname}/${dirname}_human_free.fastq
rm ${output_folder}/${dirname}/${dirname}_human_free_1.fastq
rm ${output_folder}/${dirname}/${dirname}_human_free_2.fastq
gzip ${output_folder}/${dirname}/${dirname}_human_free.fastq
# run prokka
prokka --outdir ${output_folder}/${dirname}/${dirname}_prokka_out --force --addgenes --metagenome --cpus ${thread_number} --mincontiglen 1 ${output_folder}/${dirname}/${dirname}_mh_out/final.contigs.fa
# delete megahit output folder.
#rm -r ${output_folder}/${dirname}/${dirname}_mh_out
tar -zcvf ${output_folder}/${dirname}/${dirname}_prokka_out.tar.gz ${output_folder}/${dirname}/${dirname}_prokka_out  
rm -r ${output_folder}/${dirname}/${dirname}_prokka_out
