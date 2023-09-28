#!/bin/bash

fastq-dump --split-files $1

NAME=${1##*/}
OUTFILE=${NAME%%.*}
F1=${NAME%%.*}_1.fastq
F2=${NAME%%.*}_2.fastq
CORES=$2
MEMORY_BYTES=$3

megahit -1 $F1 -2 $F2 -o ${OUTFILE}_mh_out -m $MEMORY_BYTES -t $CORES
prokka --outdir ${OUTFILE}_prokka_out --metagenome --cpus $CORES --mincontiglen 1 ${OUTFILE}_mh_out/final.contigs.fa

tar -zcvf ${OUTFILE}_prokka_out.tar.gz ${OUTFILE}_prokka_out
rm -r ${OUTFILE}_mh_out $F1 $F2 ${OUTFILE}_prokka_out
