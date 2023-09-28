#!/bin/bash

input_fastq_path=${1}
full_output=${2}
thread_number=${3}
tmp_dir=${4}

mkdir -p ${tmp_dir}
input_fastq_nopath=$(echo ${input_fastq_path##*/})
sampleName=$(echo ${input_fastq_nopath%.fastq.gz})

source activate /home/sez10/miniconda3_2/envs/metaphlan_v3.0_Humann_v3.0.0.alpha.3_curatedMetagenomicData_3.0.0_compatible

mv ${input_fastq_path} ${tmp_dir}
echo ${tmp_dir}
#rm ${input_fastq_path}

#humann --input ${input_fastq_path} --output ${full_output} --threads ${thread_number} --metaphlan-options "--index mpa_v30_CHOCOPhlAn_201901"
humann --input ${tmp_dir}/${input_fastq_nopath} --output ${tmp_dir} --threads ${thread_number} --metaphlan-options "--index mpa_v30_CHOCOPhlAn_201901"
ls -lh ${tmp_dir}
ls -lh ${tmp_dir}/${sampleName}_humann_temp
du -sh ${tmp_dir}/${sampleName}_humann_temp
du -sh ${tmp_dir}
#mv ${full_output}/${sampleName}_humann_temp/${sampleName}_metaphlan_bugs_list.tsv ${full_output}
mv ${tmp_dir}/${sampleName}_humann_temp/${sampleName}_metaphlan_bugs_list.tsv ${full_output}
mv ${tmp_dir}/${sampleName}_genefamilies.tsv ${full_output}
mv ${tmp_dir}/${sampleName}_pathabundance.tsv ${full_output}
mv ${tmp_dir}/${sampleName}_pathcoverage.tsv ${full_output}

#rm -r ${full_output}/${sampleName}_humann_temp
rm -r ${tmp_dir}

#rm ${input_fastq_path}
