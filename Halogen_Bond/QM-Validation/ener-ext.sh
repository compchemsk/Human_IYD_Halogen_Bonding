for np in ARG93 ARG97 ASN171 ASN362 ASN89 FRA440 GLU363 GLY173 GLY95 HIS348 HIS359 HIS96 LEU188 LYS356 LYS357 LYS92 MET94 THR168 TRP98 TYR360 TYR361 VAL347 VAL358 
do
rm -rf ${np}-2PT-E.dat
awk '{o = index($0, "O"); c = index($0, "C"); i = index($0, "I"); if ($0 ~ /LP/ && o && c && i && o < c && c < i) print}' ${np}-nbo.log >> ${np}-2PT-E.dat
awk '{n = index($0, "N"); c = index($0, "C"); i = index($0, "I"); if ($0 ~ /LP/ && n && c && i && n < c && c < i) print}' ${np}-nbo.log >> ${np}-2PT-E.dat
awk '{s = index($0, "S"); c = index($0, "C"); i = index($0, "I"); if ($0 ~ /LP/ && s && c && i && s < c && c < i) print}' ${np}-nbo.log >> ${np}-2PT-E.dat
cat ${np}-2PT-E.dat
done
