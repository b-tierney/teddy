#!/bin/bash

input_fastq_path=${1}
output_folder=${2}
thread_number=${3}

input_fastq_nopath=$(echo ${input_fastq_path##*/})
sampleName=$(echo ${input_fastq_nopath%_human_free.fastq.gz})
echo ${sampleName}

full_output=${output_folder}/${sampleName}
mkdir -p ${full_output}

jobName_all=$(sbatch -p transfer -t 0-0:00:05 --mem=1G /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/humann_analysis/transfer_from_standby_to_local.bash ${input_fastq_path} ${full_output})
jobName=$(echo ${jobName_all} | awk '{print $4}')

sbatch -p medium -t 1-00:00 --mem=30G -c 4 --dependency=afterok:${jobName} /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/humann_analysis/run_humann.bash ${full_output}/${input_fastq_nopath} ${full_output} ${thread_number} /tmp/${sampleName}
