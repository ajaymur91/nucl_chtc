#!/bin/bash
PROC=1
TRIALS=${2:-100000}
WDIR=$(pwd)
GRO=$1
shape=${4:-cubic}
export PLUMED_MAXBACKUP=-1
export GMX_MAXBACKUP=-1
#echo "$shape box"

mkdir -p tempdir
cd tempdir

# Place your code here 
######################################################################
# Make files

# Create Na.gro
cat > Na.gro << EOF
Na+
 1
    1Na+     Na    1   0.000   0.000   0.000
 0.0 0.0 0.0
EOF

# Create Cl.gro
cat > Cl.gro << EOF
Cl-
 1
    1Cl-     Cl    1   0.000   0.000   0.000
 0.0 0.0 0.0
EOF

# Make vdwradii.dat to hack gmx sasa
cat > vdwradii.dat << EOF
???  Na1   0.35
???  Cl1   0.35
???  Na    0.35
???  Cl    0.35
EOF

cat > contact.py << EOF
import networkx as nx
import sys
#import pylab
import copy
import pydot

# Read the contact matrix as a graph (This includes ++ and -- contacts too)
G=nx.Graph(nx.drawing.nx_pydot.read_dot(str(sys.argv[1])))

# Renumber starting from 1 instead of 0
#G=nx.relabel.convert_node_labels_to_integers(Plm, first_label=1, ordering='default', label_attribute=None)

# Remove anion-anion and cation-cation contacts from contact list
O=copy.deepcopy(G)
E=copy.deepcopy(G)
O.remove_nodes_from(tuple(str(x) for x in range(0,O.number_of_nodes()+1,2)))
E.remove_nodes_from(tuple(str(x) for x in range(1,E.number_of_nodes()+1,2)))
G.remove_edges_from(O.edges())
G.remove_edges_from(E.edges())
print(int(nx.algorithms.components.is_connected(G)))
EOF
######################################################################

cp $WDIR/$GRO $WDIR/tempdir/$GRO

# cluster -> center -> extract ions
gmx trjconv -f $GRO -o ion.gro -s $GRO -vel no &> /dev/null << EOF 
Ion
EOF

gmx trjconv -f ion.gro -s ion.gro -pbc cluster -center -o box.gro &> /dev/null << EOF
Ion
Ion
Ion
EOF

# Pad box (for inserting a pair of ions)
gmx editconf -f box.gro -bt $shape -d 0.7 -o pad.gro &> /dev/null

# insert na and cl 
gmx insert-molecules -f pad.gro -ci Na.gro -radius 0 -scale 0 -o insert1.gro -nmol 1 &> /dev/null
gmx insert-molecules -f insert1.gro -ci Cl.gro -radius 0 -scale 0 -o insert2.gro -nmol 1 &> /dev/null
echo q | gmx make_ndx -f insert2.gro &> /dev/null 
rm -rf insert1.gro insert2.gro

mkdir -p gro
echo 0 > /dev/shm/foo
for i in `seq 0 $(($TRIALS-1))`; do
	for k in `seq 1 $PROC`; do
		{ gmx insert-molecules -f pad.gro -ci Na.gro -scale 0 -o gro/"$(($i*$PROC +$k))"-1.gro -nmol 1 &> /dev/null
	       	gmx insert-molecules -f gro/"$(($i*$PROC +$k))"-1.gro -ci Cl.gro -scale 0 -o gro/"$(($i*$PROC +$k))".gro -nmol 1 &> /dev/null
		echo "
                Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
                WHOLEMOLECULES ENTITY0=Ion
                mat: CONTACT_MATRIX ATOMS=Ion SWITCH={RATIONAL R_0=0.35 NN=10000} 
                dfs: DFSCLUSTERING MATRIX=mat
                nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
                PRINT ARG=nat FILE=NAT$k
                DUMPGRAPH MATRIX=mat FILE=contact$(($i*$PROC +$k)).dot
		" | plumed driver --igro gro/$(($i*$PROC +$k)).gro --plumed /dev/stdin &> /dev/null
echo "$(($i*$PROC +$k))	  $(python contact.py contact$(($i*$PROC +$k)).dot)" >> connected
		if (( $(python contact.py contact$(($i*$PROC +$k)).dot) ))
		then 
			cp --no-clobber gro/"$(($i*$PROC +$k))".gro $WDIR/child_$GRO &> /dev/null || :
			echo 1 > /dev/shm/foo
		fi
	} &
	done
	wait
	if (( $(</dev/shm/foo) ))
	then
		break 2
	fi
done


cd $WDIR
rm -rf tempdir
######################################################################
