for i in {1..100}
do
cat Pull-${i}/Pull_dist_${i}.dat | tail -1 | awk '{print$4}' >> work-dist-kxxx.dat
done
