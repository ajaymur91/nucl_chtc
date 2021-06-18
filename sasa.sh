#!/bin/bash
PROC=1
TRIALS=${2:-500}
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
O.remove_nodes_from(tuple(str(x) for x in range(0,O.number_of_nodes()+3,2)))
E.remove_nodes_from(tuple(str(x) for x in range(1,E.number_of_nodes()+3,2)))
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
gmx editconf -f box.gro -bt $shape -d 0.35 -o pad.gro &> /dev/null

############################## SASA NA

echo 'Cl' | gmx sasa -f pad.gro -s pad.gro -probe 0 -ndots 1024 -tv V_parent_Na.xvg &> /dev/null
tail -n 1 V_parent_Na.xvg | awk '{print $2}' >> V_parent_Na
V_parent=$(tail -n 1 V_parent_Na.xvg | awk '{print $2}')

# insert Na
gmx insert-molecules -f pad.gro -ci Na.gro -radius 0 -scale 0 -o insert1.gro -nmol 1 &> /dev/null
echo q | gmx make_ndx -f insert1.gro -o index_Na.ndx &> /dev/null 
rm -rf insert1.gro 

mkdir -p gro_Na
mkdir -p child_Na
cp vdwradii.dat child_Na/

for i in `seq 0 $(($TRIALS-1))`; do
	for k in `seq 1 $PROC`; do
		{ gmx insert-molecules -f pad.gro -ci Na.gro -scale 0 -o gro_Na/"$(($i*$PROC +$k))".gro -nmol 1 &> /dev/null
		echo "$(cat << EOF
                Ion: GROUP NDX_FILE=index_Na.ndx NDX_GROUP=Ion
                WHOLEMOLECULES ENTITY0=Ion
                mat: CONTACT_MATRIX ATOMS=Ion SWITCH={RATIONAL R_0=0.35 NN=10000} 
                dfs: DFSCLUSTERING MATRIX=mat
                nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
                PRINT ARG=nat FILE=NAT
                DUMPGRAPH MATRIX=mat FILE=contact_Na$(($i*$PROC +$k)).dot
EOF
)" | plumed driver --igro gro_Na/$(($i*$PROC +$k)).gro --plumed /dev/stdin &> /dev/null
		if (( $(python contact.py contact_Na$(($i*$PROC +$k)).dot) ))
		then 
			cp gro_Na/"$(($i*$PROC +$k))".gro child_Na/"$(($i*$PROC +$k))".gro
		fi
	} &
	done
	wait
done

cd child_Na
for i in *.gro
do
	gmx editconf -f $i -bt $shape -d 0.35 -o "$i"_pad.gro &> /dev/null 
	echo "Na" | gmx sasa -f "$i"_pad.gro -s "$i"_pad.gro -probe 0 -ndots 1024 -tv volume"$i".xvg &> /dev/null
        tail -n 1 volume"$i".xvg | awk '{print $2}' >> V_child_Na
        	
done
V_child=$(awk '{ total += $1 } END { print total/NR }' V_child_Na)
echo "$V_parent * $V_child" | bc > volume_Na
#awk -v vp=$V_parent '{print $1*vp}' V_child_Na > $WDIR/V_Na
cd ..

################# SASA CL

echo 'Na' | gmx sasa -f pad.gro -s pad.gro -probe 0 -ndots 1024 -tv V_parent_Cl.xvg &> /dev/null
tail -n 1 V_parent_Cl.xvg | awk '{print $2}' >> V_parent_Cl
V_parent=$(tail -n 1 V_parent_Cl.xvg | awk '{print $2}')

# insert Cl
gmx insert-molecules -f pad.gro -ci Cl.gro -radius 0 -scale 0 -o insert1.gro -nmol 1 &> /dev/null
echo q | gmx make_ndx -f insert1.gro -o index_Cl.ndx &> /dev/null 
n=$(grep "Na\|Cl" insert1.gro | wc -l)
rm -rf insert1.gro 

mkdir -p gro_Cl
mkdir -p child_Cl
cp vdwradii.dat child_Cl/

for i in `seq 0 $(($TRIALS-1))`; do
	for k in `seq 1 $PROC`; do
		{ gmx insert-molecules -f pad.gro -ci Cl.gro -scale 0 -o gro_Cl/"$(($i*$PROC +$k))".gro -nmol 1 &> /dev/null
		echo "$(cat << EOF
                Ion: GROUP NDX_FILE=index_Cl.ndx NDX_GROUP=Ion
                WHOLEMOLECULES ENTITY0=Ion
                mat: CONTACT_MATRIX ATOMS=Ion SWITCH={RATIONAL R_0=0.35 NN=10000} 
                dfs: DFSCLUSTERING MATRIX=mat
                nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
                PRINT ARG=nat FILE=NAT
                DUMPGRAPH MATRIX=mat FILE=contact_Cl$(($i*$PROC +$k)).dot
EOF
)" | plumed driver --igro gro_Cl/$(($i*$PROC +$k)).gro --plumed /dev/stdin &> /dev/null
		sed -i "s/$(($n-1))/$n/g" contact_Cl$(($i*$PROC +$k)).dot
		if (( $(python contact.py contact_Cl$(($i*$PROC +$k)).dot) ))
		then 
			cp gro_Cl/"$(($i*$PROC +$k))".gro child_Cl/"$(($i*$PROC +$k))".gro
		fi
	} &
	done
	wait
done

cd child_Cl
for i in *.gro
do
	gmx editconf -f $i -bt $shape -d 0.35 -o "$i"_pad.gro &> /dev/null 
	echo "Cl" | gmx sasa -f "$i"_pad.gro -s "$i"_pad.gro -probe 0 -ndots 1024 -tv volume"$i".xvg &> /dev/null
        tail -n 1 volume"$i".xvg | awk '{print $2}' >> V_child_Cl
        	
done
V_child=$(awk '{ total += $1 } END { print total/NR }' V_child_Cl)
echo "$V_parent * $V_child" | bc > volume_Cl
#awk -v vp=$V_parent '{print $1*vp}' V_child_Cl > $WDIR/V_Cl
cd ..

paste child_*/volume_* | awk '{print ($1+$2)*0.5}' | tee -a volume.txt
######################################################################
# copy results back
cd $WDIR
cp $WDIR/tempdir/volume.txt $WDIR/volume_$GRO.txt
rm -rf tempdir
