#!/bin/bash

n=100

# -----------------------------------
# Step 1: Process each trajectory
# -----------------------------------
for i in $(seq 1 $n)
do
    cd Pull-${i} || exit

    rm -f rms-${i}.dat

    # Convert to DCD
    cat > cp1.in << EOF
parm hIYD-I-Tyr.parm7
trajin pull-${i}.nc
trajout pull-${i}.dcd dcd
run
quit
EOF

    cpptraj -i cp1.in

    # Extract reference structure (first frame)
    cat > cp2.in << EOF
parm hIYD-I-Tyr.parm7
trajin pull-${i}.nc 1 1 1
trajout ref.pdb pdb
run
quit
EOF

    cpptraj -i cp2.in

    rm cp2.in

    # Compute RMSD
    wordom_0.22-rc2.x86-64 -ia rmsd \
        --TITLE lig \
        --SELE "/*/@(442)/!(H*)" \
        --FIT "/*/@(1-220|221-440)/CA" \
        -imol ref.pdb \
        -itrj pull-${i}.dcd > rmsd-${i}.dat

    # Extract RMSD column
    awk '{print $2}' rmsd-${i}.dat | sed '1d' > rms-${i}.dat

    mv rms-${i}.dat ../

    # Cleanup
    rm rmsd-${i}.dat cp1.in pull-${i}.dcd ref.pdb

    cd ..
done

# -----------------------------------
# Step 2: Combine all RMSD files
# -----------------------------------
paste $(seq -f "rms-%g.dat" 1 $n) | cat -n > rmsd_full.dat

# -----------------------------------
# Step 3: Cleanup
# -----------------------------------
rm $(seq -f "rms-%g.dat" 1 $n)