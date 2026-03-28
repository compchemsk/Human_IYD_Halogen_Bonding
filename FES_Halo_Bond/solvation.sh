rm -rf solvation.dat
for np in {3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.8,10.8,11.8,12.8,13.8,14.8,15.8,16.8,17.8,18.8}
do
echo "
parm ../US-${np}/hIYD-I-Tyr.parm7
trajin ../US-${np}/${np}-us.nc 1 500 1
watershell :442 out water-${np}.dat
run
quit " > 1.in
cpptraj 1.in
awk '{sum += $2; sumsq += $2*$2; n++}
END {
    mean = sum/n
    std = sqrt(sumsq/n - mean*mean)
    print  mean, std
}' water-${np}.dat >> solvation1.dat
rm water-${np}.dat
echo "${np}" >> name.dat
done
paste name.dat solvation1.dat > solvation.dat
rm name.dat solvation1.dat
