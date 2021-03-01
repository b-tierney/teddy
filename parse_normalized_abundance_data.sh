#!/bin/bash

#for each file, split it into component parts and store those in their own little folder
for file in gene_names_batched; do

	root=$(echo $file | cut -f1 -d_)

	mkdir "${root}"_batched

	split -l 50000 $file "${root}"_batched/"${root}"_normalized_batched_

done

#create a master list of all split file locations
find . -type f -name '*batched*' > all_batch_locs

#the gene names need to be in a separate config file for the python script, so deal with them separately
grep gene_names all_batch_locs > gene_name_locs
grep -v gene_names all_batch_locs > foo
mv foo all_batch_locs

cat all_batch_locs | cut -f5 -d_ | sort | uniq > all_suffixes

#separate the master list into sublists for each suffix subset
while read line; do

	grep _"${line}"$ all_batch_locs > "${line}"_locs

done<all_suffixes


