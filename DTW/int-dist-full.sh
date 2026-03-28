#!/bin/bash

n=100

# -----------------------------
# Step 1: Extract coordinates
# -----------------------------
for i in $(seq 1 $n)
do
    cd Pull-${i} || exit

    rm -f rand_${i}.dat

    cat > cp.in << EOF
parm hIYD-I-Tyr.parm7
trajin pull-${i}.nc
vector v0 center @7314,7308,7311,7313,7309,7310 out rand_pt_${i}.dat
run
quit
EOF

    cpptraj -i cp.in

    # Extract x y z
    awk '{print $2, $3, $4}' rand_pt_${i}.dat | sed '1d' > rand_${i}.dat

    mv rand_${i}.dat ../
    rm rand_pt_${i}.dat cp.in

    cd ..
done

# -----------------------------
# Step 2: Combine all trajectories
# -----------------------------
paste $(seq -f "rand_%g.dat" 1 $n) > all_rand.dat

# -----------------------------
# Step 3: Project onto coordinate (generalized)
# -----------------------------
awk -v n=$n '
{
    for(i=1; i<=n; i++){
        x = $(3*i-2)
        y = $(3*i-1)
        z = $(3*i)

        proj = ((-58.47511066260806*x) + (75.34723814588003*y) + (52.832185591279995*z) - 3408.3603833644875) / 109.03111802401949

        printf "%f ", proj
    }
    printf "\n"
}
' all_rand.dat | cat -n > int_dis_full.dat

# -----------------------------
# Step 4: Cleanup
# -----------------------------
rm $(seq -f "rand_%g.dat" 1 $n) all_rand.dat