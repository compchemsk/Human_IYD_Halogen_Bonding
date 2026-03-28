#!/bin/bash
#install molprobity 

NCPU=40
NP_LIST="3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8"
#NP_LIST="3.8 4.0"

rm -f steric-info.dat
rm -rf tmp_parallel
mkdir tmp_parallel

#########################################
run_job() {

np=$1
n=$2

cat > tmp_parallel/tmp_${np}_${n}.in <<EOF
parm ../US-${np}/hIYD-I-Tyr.parm7
trajin ../US-${np}/${np}-us-10ns.nc ${n} ${n} 1
reference ../US-${np}/${np}-us-10ns.nc 1
strip !(:442<:5.0)
trajout tmp_parallel/${np}_${n}.pdb pdb
run
quit
EOF

cpptraj tmp_parallel/tmp_${np}_${n}.in > /dev/null 2>&1

phenix.molprobity tmp_parallel/${np}_${n}.pdb 2>/dev/null \
| grep 'All-atom Clashscore' \
| awk '{print $4}' \
> tmp_parallel/${np}_${n}.out

rm -f tmp_parallel/${np}_${n}.pdb tmp_parallel/tmp_${np}_${n}.in
}
#########################################

jobcount=0

for np in $NP_LIST
do
  for n in $(seq 1 500)
  do
    run_job $np $n &

    ((jobcount++))

    if [ "$jobcount" -ge "$NCPU" ]; then
      wait -n
      ((jobcount--))
    fi
  done
done

wait

#########################################
# Post-processing safely
#########################################

for np in $NP_LIST
do
  files=$(ls tmp_parallel/${np}_*.out 2>/dev/null)

  if [ -z "$files" ]; then
    echo "No data for $np — skipping"
    continue
  fi

  cat $files > tmp_parallel/${np}.dat

  awk '{
        sum+=$1;
        sumsq+=$1*$1
       }
       END{
        if (NR>1) {
          mean=sum/NR;
          std=sqrt((sumsq - sum*sum/NR)/(NR-1));
          printf "%.6f %.6f\n", mean, std
        }
       }' tmp_parallel/${np}.dat > tmp_parallel/${np}-stat.dat

  echo "$np" > tmp_parallel/tmp_np.dat
  paste tmp_parallel/tmp_np.dat tmp_parallel/${np}-stat.dat >> steric-info.dat

 # rm -f tmp_parallel/${np}.dat tmp_parallel/${np}-stat.dat tmp_parallel/tmp_np.dat
done

#rm -rf tmp_parallel
rm -f molprobity_coot.py molprobity.out molprobity_probe.txt

echo "Finished. Results written to steric-info.dat"
