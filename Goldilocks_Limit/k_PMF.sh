#!/bin/bash

n=100

# Step 1: Extract columns (2 and 4) from each Pull directory
for i in $(seq 1 $n)
do
    cd Pull-${i} || exit
    awk '{print $2, $4}' Pull_dist_${i}.dat > ../${i}.dat
    cd ..
done

# Step 2: Combine all files
paste $(seq -f "%g.dat" 1 $n) > data.dat

# Step 3: Compute mean (distance, force)
awk -v n=$n '{
    sum1=0; sum2=0;
    for(i=1;i<=2*n;i+=2) sum1 += $i;   # distance
    for(i=2;i<=2*n;i+=2) sum2 += $i;   # force
    print sum1/n, sum2/n
}' data.dat > mean.dat

# Step 4: Compute standard deviation (force only)
awk -v n=$n '{
    sum=0; sumsq=0;
    for(i=2;i<=2*n;i+=2){
        sum += $i;
        sumsq += $i*$i;
    }
    mean = sum/n;
    std = sqrt((sumsq/n) - (mean*mean));
    print std
}' data.dat > std.dat

# Step 5: Combine mean + std
paste mean.dat std.dat > data_PMF.dat

# Step 6: Cleanup
rm $(seq -f "%g.dat" 1 $n) data.dat mean.dat std.dat