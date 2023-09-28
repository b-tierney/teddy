#!/bin/bash

fileName=${1}

while read line
do
/n/data1/joslin/icrb/kostic/szimmerman/TEDDY_analysis/scripts/run_ttests_2.bash ${line}
done < ${fileName}
