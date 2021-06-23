#!/bin/bash
#sleep 30
# set -e
  ENVNAME=nucl
  ENVDIR=$ENVNAME
  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# these lines handle setting up the environment; you shouldn't have to modify them
  export PATH
  mkdir $ENVDIR
  tar -xzf $ENVNAME.tar.gz -C $ENVDIR
  . $ENVDIR/bin/activate
  export PLUMED_KERNEL=$SCRIPT_DIR/$ENVDIR/lib/libplumedKernel.so


NPROC=1
n=$1

##Calculate volume contribution and insert ion
 bash sasa.sh "$n".gro
 mv volume_"$n".gro.txt volume_"$n".txt
 bash insert.sh "$n".gro
 mkdir -p child
 cp child_"$n".gro child/"$n".gro

cd child
bash ../run_ti.sh $n
cd ..

# if file does not exist or empty make a files with NA
[ -s volume_"$n".txt ] || echo NA > volume_"$n".txt
[ -s FE_"$n".gro ] || echo NA > FE_"$n".gro
[ -s FE_"$n".txt ] || echo "NA NA" > FE_"$n".txt
[ -s mratio"$n" ] || echo NA > mratio"$n"

exit 0
