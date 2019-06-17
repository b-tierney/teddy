#!/bin/bash

###all_annotation_data is the catted together prokka tsv files

while read p; do

grep -w $p all_annotation_data | rev | cut -f1 | rev >> output_$1

done<$1
