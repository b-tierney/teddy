#!/bin/bash

#for each file, split it into component parts and store those in their own little folder
while read file; do root=$(echo $file | cut -f1 -d_); echo $root; awk -vRS="\n" -vORS="\t" '1'  $file > "${root}"_transposed;done<testing






awk 'BEGIN{FS="\t"; m=NUMBER }
     { for(i=1;i<=NF;++i) { 
          s = (i%m==1 ? $i : s FS $i);                                                                                                                                                 
          if (i%m==0 || i==NF) {print s > (sprintf("out.%0.5d",int(i/m)+(i%m!=0)))}
     }}' input_file