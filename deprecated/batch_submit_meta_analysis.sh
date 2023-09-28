while read p; do

a=$( echo $p | cut -f1 -d' ')
b=$( echo $p | cut -f2 -d' ')
c=$( echo $p | cut -f3 -d' ')
d=$( echo $p | cut -f4 -d' ')

sbatch -p short -c 1 --mem=50G -t 0-11:59 -n 1 ./run_meta.sh $a $b $c $d

done<todo
