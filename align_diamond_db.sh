#!/bin/bash
set -e
source activate /home/sez10/miniconda3/envs/meta_assemblers

filename=${1}
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

module load sratoolkit/2.10.7
while read dirname
do
####DOWNLOAD DIABIMMUNE FILE
wget 'https://pubs.broadinstitute.org/diabimmune/data/20/"${dirname}"_R1.fastq.gz'
wget 'https://pubs.broadinstitute.org/diabimmune/data/20/"${dirname}"_R2.fastq.gz'

####UNTAR IT APPROPRIATELY

#prefetch --force no --ngc /n/scratch3/users/b/btt5/TEDDY/prj_20109.ngc -O ${output_folder}/${dirname} $dirname
#fastq-dump --split-files --ngc /n/scratch3/users/b/btt5/TEDDY/prj_20109.ngc --outdir ${output_folder}/${dirname} ${output_folder}/${dirname}/${dirname}_dbGaP-20109.sra
# remove sra file

f1="${dirname}"_R1.fastq.gz
f2="${dirname}"_R2.fastq.gz

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
cat ${outputname} | samtools view -T ${fastafile} -@ ${thread_number} -b -h -o ${bam_filename} -
rm ${outputname}
rm ${output_folder}/${dirname}/${dirname}_human_free.fastq
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
done < ${filename}
