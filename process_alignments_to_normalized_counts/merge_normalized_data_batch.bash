#!/bin/bash
while read line
do
echo ${line}
python merge_normalized_data.py ${line}
done < ${1}
