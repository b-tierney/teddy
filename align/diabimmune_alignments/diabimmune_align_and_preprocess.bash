#!/bin/bash
set -e
source activate /home/sez10/miniconda3/envs/meta_assemblers

myfile=${1}
output_folder=${2}
adapters_reference=${3}
bitmask_file=${4}
sprism_file=${5}
temp_dir=${6}
thread_number=${7}
index=${8}
fastafile=${9}
mkdir -p ${output_folder}/${dirname}
#pairedEnd=0

export OMP_NUM_THREADS=${thread_number}

while read myline
do
dirname=$(echo "${myline}" | awk '{print $1}')
fastq_files=$(echo "${myline}" | awk '{print $2}')
fastq_one_file=$(echo ${fastq_files} | awk -F ',' '{print $1}')
fastq_two_file=$(echo ${fastq_files} | awk -F ',' '{print $2}')

echo ${dirname}
echo ${fastq_one_file}
echo ${fastq_two_file}

wget -P ${output_folder}/${dirname} ${fastq_one_file}
wget -P ${output_folder}/${dirname} ${fastq_two_file}

fastq_one_name=$(basename ${fastq_one_file})
fastq_two_name=$(basename ${fastq_two_file})

echo ${fastq_one_name}
echo ${fastq_two_name}

f1=${output_folder}/${dirname}/${fastq_one_name}
f2=${output_folder}/${dirname}/${fastq_two_name}

### check is fastq files are paired end. If not paired end script will stop
reformat.sh in1=${f1_orig} in2=${f2_orig} vpair

echo ${f1}
echo ${f2}

bbduk.sh -Xmx8g in1=${f1} in2=${f2} out1=${output_folder}/${dirname}/${dirname}_1.trimmed.fastq out2=${output_folder}/${dirname}/${dirname}_2.trimmed.fastq ref=${adapters_reference} threads=${thread_number} qin=33 ktrim=r k=23 mink=11 hdist=1 tpe tbo
# remove raw reads
rm ${f1}
rm ${f2}
# remove contamination
bmtagger.sh -b ${bitmask_file} -x ${sprism_file} -T ${temp_dir} -q1 -1 ${output_folder}/${dirname}/${dirname}_1.trimmed.fastq -2 ${output_folder}/${dirname}/${dirname}_2.trimmed.fastq -o ${output_folder}/${dirname}/${dirname}_human_free -X
rm ${output_folder}/${dirname}/${dirname}_1.trimmed.fastq
rm ${output_folder}/${dirname}/${dirname}_2.trimmed.fastq

## cat the 2 files together for aligning to diamond
cat ${output_folder}/${dirname}/${dirname}_human_free_1.fastq ${output_folder}/${dirname}/${dirname}_human_free_2.fastq > ${output_folder}/${dirname}/${dirname}_human_free.fastq
rm ${output_folder}/${dirname}/${dirname}_human_free_1.fastq
rm ${output_folder}/${dirname}/${dirname}_human_free_2.fastq

bam_filename=${output_folder}/${dirname}/${dirname}'.catalog.bam'
outputname=${output_folder}/${dirname}/${dirname}_output
echo "Starting Alignment." >&2
diamond blastx --db ${index} --query ${output_folder}/${dirname}/${dirname}_human_free.fastq --outfmt 101 -p ${thread_number} -o ${outputname}
gzip ${output_folder}/${dirname}/${dirname}_human_free.fastq
cat ${outputname} | samtools view -T ${fastafile} -@ ${thread_number} -b -h -o ${bam_filename} -
rm ${outputname}
echo 'Sorting the bam file' >&2
samtools 'sort' \
    -l 9 \
    -o ${bam_filename%.*}'.sorted.bam' \
    -O bam \
    -@ ${thread_number} \
    ${bam_filename}

# Cleaning up unsorted bam
echo 'Removing unsorted bam' >&2
rm ${bam_filename}
bam_filename=${bam_filename%.*}'.sorted.bam'

# Index the bam
echo 'Indexing the bam file' >&2
bam_index_filename=${bam_filename%.*}'.bai'
samtools 'index' -@ ${thread_number} -b ${bam_filename} ${bam_index_filename}
echo 'Extracting aligned reads'
countsfile=${output_folder}/${dirname}/${dirname}_alignment_data.tsv
samtools idxstats -@ ${thread_number} ${bam_filename} > ${countsfile}
# remove last line of output file cause it has number of unaligned reads
head -n -1 ${countsfile} > ${output_folder}/${dirname}/${dirname}_temporary
mv ${output_folder}/${dirname}/${dirname}_temporary ${countsfile}

rm ${bam_index_filename}
rm ${bam_filename}
gzip ${countsfile}
done < ${myfile}
