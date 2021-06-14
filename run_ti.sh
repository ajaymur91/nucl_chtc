#!/bin/bash


#find . -name "#*#" -exec rm -rf {} \;
NSTEPS=100000 # NSTEPS/2000 ps
NPROC=1
GRO="$1".gro
WDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

#########################################################################
# Place your code here 
######################################################################

################### Make itp ################
mkdir -p itp
cat << EOF > itp/template.top
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5           0.833333

[ atomtypes ]
; name    at.num    mass    charge ptype  sigma      epsilon
Na1           11  22.989769  0.00000000  A         0.2584         0.4184
Cl1           17  35.453200  0.00000000  A         0.4036         0.4184
O1             8  15.999430  0.00000000  A     0.31657195      0.6497752
H1             1   1.007947  0.00000000  A          0.065        0.16628


[ moleculetype ]
; Name            nrexcl
Na+          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Na+ rtp Na+ q 1.0
    1        Na1      1    Na+     Na      1 0.50000000  22.989769   ; qtot 1.000000


[ moleculetype ]
; Name            nrexcl
Cl-          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Cl- rtp Cl- q -1.0
    1        Cl1      1    Cl-     Cl      1 -0.50000000  35.453200   ; qtot -1.000000

[ moleculetype ]
; Name            nrexcl
ghost          0

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Na+ rtp Na+ q 1.0
    1        Na1      1    Na+     Na      1 0.50000000  22.989769   ; qtot 1.000000
    2        Cl1      2    Cl-     Cl      2 -0.50000000  35.453200   ; qtot -1.000000
[ moleculetype ]
; Name            nrexcl
HOH          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 HOH rtp HOH q 0.0
    1         O1      1    HOH      O      1 -0.84760000  15.999430   ; qtot -0.847600
    2         H1      1    HOH     H1      2 0.42380000   1.007947   ; qtot -0.423800
    3         H1      1    HOH     H2      3 0.42380000   1.007947   ; qtot 0.000000

#ifdef FLEXIBLE

[ bonds ]
;    ai     aj funct         c0         c1         c2         c3
      1      2     1   0.10000 462750.400000
      1      3     1   0.10000 462750.400000

[ angles ]
;    ai     aj     ak funct         c0         c1         c2         c3
      2      1      3     1   109.4700000 836.800000


#else

[ settles ]
; i     funct   doh     dhh
1     1   0.10000000   0.16329809

#endif

[ exclusions ]
1  2  3
2  1  3
3  1  2

[ system ]
; Name
Generic title

[ molecules ]
; Compound       #mols
EOF
################### Make mdp ################
mkdir -p mdp

cat << EOF > mdp/min.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1.0          ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
cutoff-scheme            = group    ; Buffered neighbor searching
rlist                    = 0   ;gas ph min (0 means no cutoff)
rcoulomb                 = 0   ;gas ph min (0 means no cutoff)
rvdw                     = 0   ;gas ph min (0 means no cutoff)
pbc                      = no
nstenergy                = 10
nstxout                  = 10
nstlist                  = 0
ns-type                  = simple
continuation             = no  ; does the same thing as unconstrained_start
;unconstrained_start     = yes ; depricated

;       Vacuum simulations are a special case, in which neighbor lists and cutoffs are basically thrown out.
;       All interactions are considered (rlist = rcoulomb = rvdw = 0) and the neighbor list is fixed (nstlist = 0).
EOF

cat << EOF > mdp/md.mdp
comm-mode                = Linear
nstcomm                  = 10

integrator              = md        ; leap-frog integrator
nsteps                  = 10000000    ; 2 * 1000000 = 20000 ps 
dt                      = 0.001     ; 2 fs

nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 0      ; save energies every 10.0 ps
nstlog                  = 10000      ; update log file every 10.0 ps
nstxout-compressed      = 1000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system

continuation            = no       ; Restarting after NPT 
;gen-vel                 = yes
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)

coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

tcoupl                  = nose-hoover           ; modified Berendsen thermostat
tc-grps                 = system                ; two coupling groups - more accurate
tau_t                   = 2.0                     ; time constant, in ps
ref_t                   = 300                   ; reference temperature, one for each group, in K

;Pressure coupling is on
;pcoupl                  = parrinello-rahman     ; Pressure coupling on in NPT
;pcoupltype              = isotropic             ; uniform scaling of box vectors
;tau_p                   = 4.0                   ; time constant, in ps
;ref_p                   = 1.0                   ; reference pressure, in bar
;compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
;refcoord_scaling        = com
pbc                     = xyz       ; 3-D PBC
; DispCorr                = EnerPres  ; account for cut-off vdW scheme

;--------------------
; Free energy parameters
free-energy               = yes

sc-alpha                  = 0.5 
sc-r-power                = 6 
sc-power                  = 1  

init-lambda-state        = XXX
coul-lambdas             = 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
vdw-lambdas              = 0.0 0.2 0.4 0.6 0.8 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0


calc-lambda-neighbors = -1
nstdhdl                  = 1 
dhdl-print-energy        = yes
couple-moltype           = ghost
couple-lambda0           = none
couple-lambda1           = vdw-q
; we are keeping the intramolecular interactions ON in all the interactions from state 0 to state 8
couple-intramol          = yes
EOF

 

cat << EOF > mdp/npt.mdp
comm-mode                = Linear
nstcomm                  = 10

integrator              = md        ; leap-frog integrator
nsteps                  = 10000000    ; 2 * 1000000 = 20000 ps 
dt                      = 0.001     ; 2 fs

nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 0      ; save energies every 10.0 ps
nstlog                  = 0      ; update log file every 10.0 ps
nstxout-compressed      = 1000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system

gen-vel                 = yes
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy

cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)

coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT

tcoupl                  = nose-hoover           ; modified Berendsen thermostat
tc-grps                 = system                ; two coupling groups - more accurate
tau_t                   = 2.0                     ; time constant, in ps
ref_t                   = 300                   ; reference temperature, one for each group, in K

pbc                     = xyz       ; 3-D PBC

;Pressure coupling is on
;pcoupl                  = parrinello-rahman     ; Pressure coupling on in NPT
;pcoupltype              = isotropic             ; uniform scaling of box vectors
;tau_p                   = 4.0                   ; time constant, in ps
;ref_p                   = 1.0                   ; reference pressure, in bar
;compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
;refcoord_scaling        = com
; DispCorr                = EnerPres  ; account for cut-off vdW scheme
EOF
###########################################################

GMX=gmx
GMX_MPI=gmx
#MPIRUN_1="mpirun -n 1"
MPIRUN_1=" "

n_ti=16

mkdir -p run_ti
N_ions=$(grep "Na" $GRO | wc -l)
N_plain=$((N_ions-1))
cp itp/template.top run_ti/ghost.top

for i in `seq 1 $N_plain`
do
cat << EOF >> run_ti/ghost.top
Na+                  1
Cl-                  1
EOF
done
cat << EOF >> run_ti/ghost.top
ghost                1
EOF
cp $GRO run_ti/start.gro


echo 'q' | $GMX make_ndx -f $GRO &> /dev/null

   echo "$(cat << EOF
        Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion
        WHOLEMOLECULES ENTITY0=Ion
        mat: CONTACT_MATRIX ATOMS=Ion SWITCH={RATIONAL R_0=0.35 NN=10000} NOPBC
        dfs: DFSCLUSTERING MATRIX=mat
        nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
        PRINT ARG=nat FILE=NAT
        DUMPGRAPH MATRIX=mat FILE=contact.dot
EOF
)" | plumed --no-mpi driver --igro $GRO --plumed /dev/stdin &> /dev/null;


#############
cat << EOF > contact.py
import networkx as nx
import sys
#import pylab
import copy
import pydot

# Read the contact matrix as a graph (This includes ++ and -- contacts too)
G=nx.Graph(nx.drawing.nx_pydot.read_dot('contact.dot'))

# Remove anion-anion and cation-cation contacts from contact list
O=copy.deepcopy(G)
E=copy.deepcopy(G)
O.remove_nodes_from(tuple(str(x) for x in range(0,O.number_of_nodes()+3,2)))
E.remove_nodes_from(tuple(str(x) for x in range(1,E.number_of_nodes()+3,2)))
G.remove_edges_from(O.edges())
G.remove_edges_from(E.edges())
print(int(nx.algorithms.components.is_connected(G)))


sys.stdout = open("mst", "w")
if nx.algorithms.components.is_connected(G):
  MST=nx.algorithms.tree.minimum_spanning_tree(G)
  for line in nx.generate_edgelist(MST, data=False):print(line)

sys.stdout.close()
EOF
#############

   if (( $(python contact.py contact.dot) ))
   	   then
           # Create plumed.dat to implement MST
           echo "Ion: GROUP NDX_FILE=index.ndx NDX_GROUP=Ion" > plumed.dat
           echo -e "WHOLEMOLECULES ENTITY0=Ion\n" >> plumed.dat

           n=1
           awk '{print $1+1","$2+1}' mst > edges
           while read p
           do
           echo "#D($p)" >> plumed.dat
           echo "DISTANCE ATOMS=$p LABEL=d$n NOPBC" >> plumed.dat
           echo "UPPER_WALLS ARG=d$n AT=0.35 KAPPA=1000.0 EXP=2 EPS=1 OFFSET=0 LABEL=uwall$n" >> plumed.dat
           echo " " >> plumed.dat
           n=$((n+1))
           done < edges
   fi
##############

# Make box bigger (does not really matter - but do it anyway to be safe)
$GMX editconf -f $GRO -bt cubic -box 4.2139 4.2139 4.2139 -o box.gro # &> /dev/null
# $GMX solvate -scale 0.4 -cp box.gro -cs struct/spc216.gro -p run_ti/ghost.top -o solv.gro -maxsol 2500 # &> /dev/null

   $MPIRUN_1 $GMX_MPI grompp -c box.gro  -o min.tpr -f mdp/min.mdp -p run_ti/ghost.top
   time $MPIRUN_1 $GMX_MPI mdrun -v -deffnm min -plumed plumed.dat

   $MPIRUN_1 $GMX_MPI grompp -c min.gro  -o npt.tpr -f mdp/npt.mdp -p run_ti/ghost.top
   time $MPIRUN_1 $GMX_MPI mdrun -v -deffnm npt -plumed plumed.dat -ntomp 1 -nsteps $NSTEPS

   cp index.ndx run_ti/
   cp plumed.dat run_ti/
   sed -i 's/index/\.\.\/index/g' run_ti/plumed.dat

for i in `seq 0 $(($n_ti-1))`; do echo $i; sed "s/XXX/$i/g" mdp/md.mdp > run_ti/sd."$i".mdp ; done
for i in `seq 0 $(($n_ti-1))`; do echo $i; mkdir -p run_ti/topol"$i"; $MPIRUN_1 $GMX_MPI grompp -f run_ti/sd."$i".mdp -c npt.gro -p run_ti/ghost.top -o run_ti/topol"$i".tpr -maxwarn 1 ; done

cd run_ti
n_ti=16
i=0
N=$NPROC
echo "N = $N, n_ti = $n_ti"
(
for j in `seq 0 $(($n_ti-1))`; do
   ((i=i%N)); ((i++==0)) && wait
   { echo $j; $GMX_MPI mdrun -ntomp 1 -deffnm topol"$j" -dhdl topol."$j".dhdl.xvg -nsteps $NSTEPS -plumed ../plumed.dat &> /dev/null; }&
done
wait
)


mkdir -p dhdl 
cp *.dhdl.xvg dhdl/
cd dhdl
alchemical_analysis -m 'ti+bar+mbar' -p topol -s $(($NSTEPS/2000)) -v
cp results.txt $WDIR/FE_"$1".txt
cd ..

echo 0 | gmx trjconv -f topol$(($n_ti-1)).xtc -s topol$(($n_ti-1)).tpr -o topol$(($n_ti-1))/md.gro -sep -b $(($NSTEPS/2000)) -dt 1
cd topol$(($n_ti-1))
ls *.gro | xargs -P 1 -I{} bash $WDIR/extract_connected.sh {}
ls is_connected*.gro | xargs -P 1 -I{} bash $WDIR/M_ratio.sh {} | awk '{ total += $1; count++ } END { print total/count }' > $WDIR/mratio"$1"
ls is_connected* | sort -V | head -n 1 | xargs -I{} cp "{}" $WDIR/FE_"$1".gro
cd ..

cd ..

cd $WDIR 
