for np in ARG93 ARG97 ASN171 ASN362 ASN89 FRA440 GLU363 GLY173 GLY95 HIS348 HIS359 HIS96 LEU188 LYS356 LYS357 LYS92 MET94 THR168 TRP98 TYR360 TYR361 VAL347 VAL358
do
echo "${np}" > rol1.dat	
cat ${np}.dat > rol2.dat
paste rol1.dat rol2.dat >> role1.dat
rm rol1.dat rol2.dat
done
paste role1.dat 2d-halo-ener.dat > role-halo.dat
rm role1.dat
