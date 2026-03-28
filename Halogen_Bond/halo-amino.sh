grep -v 'For trajectory' halo.out > 1.dat
grep -v 'Frame' 1.dat > 2.dat
grep -v 'No halogen' 2.dat > 3.dat
cat 3.dat | sort -nk1 > 4.dat
grep -v 'HOH' 4.dat > 5.dat
grep -v 'HOH' 4.dat | uniq > 5.dat
grep -v 'HOH' 4.dat  > 5.dat
awk '{print$1}' 5.dat > 6.dat
cat 6.dat | uniq > halo-amino.dat
rm 1.dat 2.dat 3.dat 4.dat 5.dat 6.dat 
