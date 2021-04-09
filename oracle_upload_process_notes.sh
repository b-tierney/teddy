oci os object bulk-upload -bn t1d_predictioni --src-dir healthy_pre-sero-3month --dry-run

https://objectstorage.eu-zurich-1.oraclecloud.com/p/_EUc6aRM5WFhjDPyZ2WtsufVBQ3vzYwQDTAddyaPN_vgZvLI1xxkIKaGgB7aBig2/n/zrjwsvatolrg/b/t1d_prediction/o/


find mydir -type f -exec curl -u xxx:psw --ftp-create-dirs -T {} ftp://192.168.1.158/public/demon_test/{} \;


for file in $1/*;

do 

echo curl -X PUT --data-binary @"$file" https://objectstorage.eu-zurich-1.oraclecloud.com/p/_EUc6aRM5WFhjDPyZ2WtsufVBQ3vzYwQDTAddyaPN_vgZvLI1xxkIKaGgB7aBig2/n/zrjwsvatolrg/b/t1d_prediction/o/"$file"

done


while read p; do 

p=$(echo $p | sed 's/.\///g')

echo $p >> output_counts

oci os object list --all -bn t1d_prediction --prefix $p | grep name | wc -l >> output_counts

done<uploadfolds