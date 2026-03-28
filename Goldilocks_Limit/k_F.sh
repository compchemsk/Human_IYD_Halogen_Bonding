#!/bin/bash

# Extract data from each Pull directory
for i in {1..100}
do
    cd Pull-${i} || exit
    awk '{print $2, $3}' Pull_dist_${i}.dat > ${i}.dat
    mv ${i}.dat ../
    cd ..
done

# Combine all 100 files
paste {1..100}.dat > data.dat

# Compute average of columns
awk '{
    sum1=0; sum2=0;
    for(i=1;i<=199;i+=2) sum1 += $i;
    for(i=2;i<=200;i+=2) sum2 += $i;
    print sum1/100, sum2/100
}' data.dat > data_force.dat

# Cleanup
rm {1..100}.dat data.dat