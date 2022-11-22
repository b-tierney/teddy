#!/bin/bash
#
set -e
#source activate /home/btt5/miniconda3
source activate /home/sez10/miniconda3/envs/meta_assemblers
num_threads=4

filename=${1}
output_dir=${2}
index=${3}
fastafile=${4}
while read name
do
basename="$(echo $name | rev | cut -f2 -d/ | rev)"
output_subdir=${output_dir}/${basename}
outputname=${output_subdir}/"${basename}"_output.gz
mkdir -p ${output_subdir}

#cat ${name}/* > "${name}""${basename}"_temp

bam_filename=${output_subdir}/${basename}'.catalog.bam'
echo "Starting Alignment." >&2
diamond blastx --db ${index} --query "${name}" --outfmt 101 -p ${num_threads} --compress 1 --unal 0 -o ${outputname}
zcat ${outputname} | samtools view -T ${fastafile} -@ ${num_threads} -b -h -o ${bam_filename} -

rm ${outputname}
#rm "${name}""${basename}"_temp
#
# Cleaning up fastq's to make space
#echo 'Removing FastQ files' >&2
#rm ${read1_filename} ${read2_filename}

#
# Sort the bam file
echo 'Sorting the bam file' >&2
samtools 'sort' \
    -l 9 \
    -o ${bam_filename%.*}'.sorted.bam' \
    -O bam \
    -@ ${num_threads} \
    ${bam_filename}

#
# Cleaning up unsorted bam
echo 'Removing unsorted bam' >&2
rm ${bam_filename}
bam_filename=${bam_filename%.*}'.sorted.bam'

#
# Index the bam
echo 'Indexing the bam file' >&2
bam_index_filename=${bam_filename%.*}'.bai'
samtools 'index' -@ ${num_threads} -b ${bam_filename} ${bam_index_filename}

echo 'Extracting aligned reads'
countsfile=${output_subdir}/${basename}_alignment_data.tsv
samtools idxstats -@ ${num_threads} ${bam_filename} > ${countsfile}
rm ${bam_filename}
rm ${bam_index_filename}
head -n -1 ${countsfile} > ${countsfile}_temporary
mv ${countsfile}_temporary ${countsfile}
gzip ${countsfile}
done < ${filename}
#scp ${countsfile} ${countsfile}_raw_counts.tsv
#echo 'Normalizing'
#awk '{printf "%.20f\t",$3 / $2}1' ${countsfile} > ${basename}_countstemp
#mv ${basename}_countstemp ${countsfile}
#totalaligncount=$(cut -f1 ${countsfile} | awk '{s+=$1} END {print s}')
#awk '{printf "%.20f\t",$1/"'${totalaligncount}'"}1' ${countsfile} > ${basename}_countstemp
#cut -f1 ${basename}_countstemp > ${countsfile}

#sed -i 's/0.00000000000000000000/0.0/g' ${countsfile}
#rm ${countsfile}
#pigz ${bam_index_filename}
#pigz ${bam_filename}
#rm -r ${basename}_countstemp
#${countsfile}_raw_counts.tsv

# Finished.
#echo 'Finished' >&2
#exit 0

