#!/bin/bash
if [[ $# -ne 1 ]]; then
        echo " Pass the configuration (gro format)" >&2
        echo " Usage: bash extract_connected.sh <config.gro>" >&2
        exit 2
fi

GRO=$1
WDIR=$(pwd)
tempdir=$(mktemp -dt "$(basename $0).XXXXXXXX" --tmpdir=$XDG_RUNTIME_DIR)

cp $GRO $tempdir/
cd $tempdir
#########################################################################
# Place your code here 
######################################################################

# Make index file
echo 'q' | gmx make_ndx -f "$GRO" &> /dev/null
echo "$(cat << EOF
        Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
        WHOLEMOLECULES ENTITY0=Ion
        mat: CONTACT_MATRIX ATOMS=Ion SWITCH={RATIONAL R_0=0.35 NN=10000} NOPBC
        dfs: DFSCLUSTERING MATRIX=mat
        nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
        PRINT ARG=nat FILE=NAT
        DUMPGRAPH MATRIX=mat FILE=contact.dot
EOF
)" | plumed driver --igro "$GRO" --plumed /dev/stdin &> /dev/null;

#####################################################################
cat << EOF > contact.py
import networkx as nx
import sys
#import pylab
import copy
#import pydot

# Read the contact matrix as a graph (This includes ++ and -- contacts too)
G=nx.Graph(nx.drawing.nx_pydot.read_dot('contact.dot'))
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
#######################################################################
rv=$(python contact.py contact.dot)
if (( $rv ))
	then
	cp "$GRO" $WDIR/is_connected_"$GRO"
fi

# copy results back
cd $WDIR
rm -rf $tempdir
echo $rv
