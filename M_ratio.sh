#!/bin/bash

GRO=$1
WDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

mkdir -p tempdir
cp $GRO tempdir/
cd tempdir
#########################################################################
# Place your code here 
######################################################################

# make index file
echo q | gmx make_ndx -f $GRO &> /dev/null

# pbc Wrap cluster 
echo "$(cat << EOF 
Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
WRAPAROUND ATOMS=Ion AROUND=1
DUMPATOMS FILE=dump.gro ATOMS=Ion
EOF
)" | plumed driver --igro $GRO --plumed /dev/stdin &> /dev/null 

# Center cluster in box
gmx trjconv -f dump.gro -s $GRO -center -o box.gro &> /dev/null << EOF
Ion
Ion
EOF

# Pad box
gmx editconf -f box.gro -bt cubic -d 0.7 -o pad.gro &> /dev/null

# Create matrix.dat - plumed creates a contact matrix for the cluster
cat > matrix.dat << EOF
Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
mat: CONTACT_MATRIX ATOMS=Ion SWITCH={RATIONAL R_0=0.35 NN=10000}
dfs: DFSCLUSTERING MATRIX=mat
nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
PRINT ARG=nat FILE=NAT
DUMPGRAPH MATRIX=mat FILE=matrix.dot
EOF

# Create a contact matrix and save to matrix.dot 
plumed driver --igro pad.gro --plumed matrix.dat &> /dev/null

# Make M_ratio.py: Python script that reads matrix.dot and calculates M_ratio
cat << EOF > M_ratio.py
import networkx as nx
import copy
import pydot
# Read the contact matrix as a graph (This includes ++ and -- contacts too)
G=nx.Graph(nx.drawing.nx_pydot.read_dot("./matrix.dot"))

# Remove anion-anion and cation-cation contacts from contact list
O=copy.deepcopy(G)
E=copy.deepcopy(G)
O.remove_nodes_from(tuple(str(x) for x in range(0,O.number_of_nodes()+3,2)))
E.remove_nodes_from(tuple(str(x) for x in range(1,E.number_of_nodes()+3,2)))
G.remove_edges_from(O.edges())
G.remove_edges_from(E.edges())

# n: Track the number of (N-1) sub-clusters satisfying Stillinger criteria
n=0
for enode in G.nodes():
    if (int(enode) % 2 == 0):
        for onode in G.nodes():
            if (int(onode) % 2 == 1):
                #print(enode,onode)
                H = copy.deepcopy(G)
                H.remove_node(enode)
                H.remove_node(onode)
                #Check to see if the (N-1) sub-cluster is connected
                if nx.algorithms.components.connected.is_connected(H): n += 1

M=(G.number_of_nodes()**2)/float(4*n)
print("{:.5f}".format(M))
EOF

# Run the python script to calculate M_ratio and print to stdout
python M_ratio.py

######################################################################
# copy results back
cd $WDIR
rm -rf tempdir

