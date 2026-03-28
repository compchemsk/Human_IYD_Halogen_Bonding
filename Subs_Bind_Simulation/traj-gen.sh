#!/bin/bash

PARM="I-Tyr-binding.parm7"
TRAJ="sim1.nc"

OFFSET=2419

for c in {0..5}
do

INDEXFILE="cluster_${c}_nearest100.dat"

echo "Processing cluster $c"

# ==========================
# Generate frame numbers
# ==========================

awk -v off=$OFFSET 'NR>1 {print $1 + off}' $INDEXFILE > frames.dat


# ==========================
# Extract frames as PDB
# ==========================

count=0

while read frame
do

count=$((count+1))

cpptraj <<EOF > /dev/null
parm $PARM
trajin $TRAJ $frame $frame
trajout frame_${count}.pdb pdb include_ep
run
quit
EOF

done < frames.dat

echo "Cluster $c PDB frames created"


# ==========================
# Merge PDBs into XTC
# ==========================

cpptraj <<EOF > /dev/null
parm $PARM
trajin first.pdb
trajin frame_*.pdb
trajout cluster${c}.xtc xtc include_ep
run
quit
EOF

echo "cluster${c}.xtc created"


# ==========================
# Clean temporary files
# ==========================

rm frame_*.pdb
rm frames.dat

done

echo "All clusters finished."
