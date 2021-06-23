#!/bin/bash
#set -e
#if [[ $# -ne 3 ]]; then
#        echo " Pass the cluster size and trials and Number of samples" >&2
#        echo " Usage: bash generate_child.sh <cluster_size> <trials> <samples>" >&2
#        exit 2
#fi
eval "$(conda shell.bash hook)"
conda activate NUCL

N=$(($1+1))
T=100
Ns=100
Gm=0.1

eval echo volume_{0..$((T-1))}.txt | xargs -n 1 cat > Veff.txt
eval echo mratio{0..$((T-1))} | xargs -n 1 cat > mratio.txt
eval echo FE_{0..$((T-1))}.txt | xargs -n 1 tail -n 1 -q | awk '{print $2}' > TI.txt

[ -s Pr.txt ] || for i in `seq 1 $Ns`; do     echo 1; done > Pr.txt
echo "Gm=$Gm"

Rscript --vanilla Boltzmann.R $N $Ns $Gm
mkdir -p N_$(($N+1))

k=0

while read p; do
#echo "$p, $k"
  	cp FE_"$p".gro N_$(($N+1))/"$k".gro
  	k=$((k+1))
done <boltzmann.txt

mkdir -p results/N$N

mv TI.txt results/N$N/
mv FE*.txt results/N$N/
mv FE*.gro results/N$N/
mv mratio.txt results/N$N/
mv mratio* results/N$N/
mv Veff.txt results/N$N/
mv volume*.txt results/N$N/
mv boltzmann.txt results/N$N/
mv dG.txt results/N$N/
cp Pr.txt results/N$N/
cp P.txt results/N$N/


RETRY=$1
if (( $RETRY < 10 ))
then {
   echo "GO AGAIN"
   exit 1
}
else {
   echo "DONE"
   exit 0
}
fi
