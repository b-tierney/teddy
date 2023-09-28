#!/bin/bash

while read p; do while read pp; do echo Rscript prep_for_voe.R $p $pp sampling_mapping_files; done<metadata_file_names; done<batched_file_names > all_jobs

split -l2 all_jobs average_

sed -i "1i #!/bin/bash" average_*

chmod +x average_*

for file in average*; do echo sbatch -p short -c 1 -n 1 --mem=40G -t 0-11:59 ./"${file}"; done > all_submission_scripts

split -l500 all_submission_scripts average_jobs_sbatch_

chmod +x average_jobs_sbatch_*

#./average_jobs_sbatch_aa 
#./average_jobs_sbatch_ab  
#./average_jobs_sbatch_ac  
#./average_jobs_sbatch_ad