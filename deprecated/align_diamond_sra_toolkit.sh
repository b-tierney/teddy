#!/bin/bash
#

date > log_${1}
####index AND azure batch run + prep code AND indexed catalog for view in the provisioning folder!!
#!/bin/bash

module load gcc
module load sratoolkit
module load samtools
module load bowtie2

prefetch $1
fastq-dump --split-files sra/$1'.sra'
num_threads=$2
#actually run the alignment
fasta_filename=$1
read1_filename=${fasta_filename}_1.fastq
read2_filename=${fasta_filename}_2.fastq

index=${3}

cat ${read1_filename} ${read2_filename} > foo_${fasta_filename}

bam_filename=${fasta_filename}'.catalog.bam'
echo "Starting Alignment." >&2
diamond blastx --db ${index} --query foo_${fasta_filename} --outfmt 101 -b 6 -p ${num_threads} -o output_${fasta_filename} 

cat output_${fasta_filename} | samtools view -T ${4} -@ ${num_threads} -b -h -o ${bam_filename} -

rm output_${fasta_filename} 
rm foo_${fasta_filename} 
rm ${read1_filename} ${read2_filename}

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
countsfile=${fasta_filename}_alignment_data.tsv
samtools idxstats -@ ${num_threads} ${bam_filename} > ${countsfile}
head -n -1 ${countsfile} > temp
mv temp ${countsfile}

scp ${countsfile} ${fasta_filename}_raw_counts.tsv
echo 'Normalizing'
awk '{printf "%.20f\t",$3 / $2}1' ${countsfile} > foo_${fasta_filename} 
mv foo_${fasta_filename}  ${countsfile}
totalaligncount=$(cut -f1 ${countsfile} | awk '{s+=$1} END {print s}')
awk '{printf "%.20f\t",$1/"'${totalaligncount}'"}1' ${countsfile} > foo_${fasta_filename} 
echo $fasta_filename > ${countsfile}
cut -f1 foo_${fasta_filename}  >> ${countsfile}

sed -i 's/0.00000000000000000000/0.0/g' ${countsfile}

rm foo_${fasta_filename} 

# Finished.
echo 'Finished' >&2
exit 0

date >> log
