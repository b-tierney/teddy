#!/bin/bash

#config for oracle

#download data
oci os object bulk-download --overwrite -bn t1d_prediction --include "${1}"/"${1}"* --download-dir .

mv run_regressions* "${1}"/

mv queue.py "${1}"/

cd "${1}"

#generate data locs and commands file

ls healthy* | grep -v metadata > datalocs

while read p; do echo Rscript run_regressions.R $p *metadata*; done<datalocs > commands

python queue.py