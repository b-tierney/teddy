oci os object bulk-upload -bn t1d_prediction/healthy_pre-t1d-24month --src-dir healthy_pre-t1d-24month --dry-run

https://objectstorage.eu-zurich-1.oraclecloud.com/p/_EUc6aRM5WFhjDPyZ2WtsufVBQ3vzYwQDTAddyaPN_vgZvLI1xxkIKaGgB7aBig2/n/zrjwsvatolrg/b/t1d_prediction/o/


find mydir -type f -exec curl -u xxx:psw --ftp-create-dirs -T {} ftp://192.168.1.158/public/demon_test/{} \;


for file in $1/*;

do 

echo curl -X PUT --data-binary @"$file" https://objectstorage.eu-zurich-1.oraclecloud.com/p/97g1B4gny3r51Yghs-sYwrXjpLOOZo8_fLqXVQeer5ycaLFF8Z6MgXlNNrbGd4kH/n/zrjwsvatolrg/b/t1d_prediction/o/"$file"

done


while read p; do 

p=$(echo $p | sed 's/.\///g')

echo $p >> output_counts

oci os object list --all -bn t1d_prediction --prefix $p | grep name | wc -l >> output_counts

done<uploadfolds


curl -X PUT --data-binary @healthy_any-T1D/full_association_output_adjusted.rds https://objectstorage.eu-zurich-1.oraclecloud.com/p/97g1B4gny3r51Yghs-sYwrXjpLOOZo8_fLqXVQeer5ycaLFF8Z6MgXlNNrbGd4kH/n/zrjwsvatolrg/b/t1d_prediction/o/healthy_any-T1D/full_association_output_adjusted.rds

oci os object list -bn t1d_prediction --limit 10000  | grep name | cut -f1 -d/ | uniq -c | sort -n


curl -X get --data-binary https://objectstorage.eu-zurich-1.oraclecloud.com/p/97g1B4gny3r51Yghs-sYwrXjpLOOZo8_fLqXVQeer5ycaLFF8Z6MgXlNNrbGd4kH/n/zrjwsvatolrg/b/t1d_prediction/o/healthy_HLA/full_association_output_adjusted.rds @healthy_HLA/full_association_output_adjusted.rds 

#!/bin/bash

#config for oracle

#download data
oci os object bulk-download --overwrite -bn t1d_prediction --include "${1}"/"${1}"* --download-dir .