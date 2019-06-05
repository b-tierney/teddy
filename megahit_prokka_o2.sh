#!/bin/bash

#megahit/prokka deployment on TEDDY data

FILE1=$1

while read p; do

sbatch -p short -t 0-11:59 -c 1 --mem=6G ./run_megahit_and_prokka_teddy.sh $p 1 60000000000

done<$FILE1
