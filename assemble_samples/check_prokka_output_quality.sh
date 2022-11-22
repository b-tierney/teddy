#!/bin/bash

filename=$(ls $1/*faa)
filename2=$(echo $1 |rev  | cut -f1 -d'/' | rev )

echo $filename
echo $filename2

rm -rf prokka_quality_data/${filename2}_data

filename3=$(echo $filename |rev  | cut -f1 -d'/' | rev )
num=$(grep '>' "${filename}"|wc -l)
typeof=$(file "${filename}" | cut -f2- -d' ')
blanks=$(grep -cvP '\S' "${filename}")
toomany=$(grep '>>' "${filename}" | wc -l)
seqid=$(head -1 "${filename}" | cut -f1 -d'_' | sed 's/>//g')
nonascii=$(grep -P -n "[\x80-\xFF]" "${filename}"|wc -l)
echo -e "${filename3}""\t""${filename2}""\t""${toomany}""\t""${typeof}""\t""${seqid}""\t""${num}""\t""${blanks}""\t""${nonascii}" >> prokka_quality_data/${filename2}_data



