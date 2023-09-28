#!/bin/bash

myfasta=${1}
output=${2}
/home/sez10/kostic_lab/software/diamond_v2.0.11/diamond blastp -d /n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/gencode.v40.pc_translations_db -q ${myfasta} -o ${output} -p 4 --ultra-sensitive
