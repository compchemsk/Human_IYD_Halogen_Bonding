#!/bin/bash
set -e
#This is an example code only. Actual code will change based on how much is the length of the simulation.

################ USER SETTINGS ################
NPWHAM=./npwham
T=300

windows=(3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 \
         6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 \
         9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8)

k_high=10.0   # kcal/mol/Å^2
k_low=2.5
################################################

nw=${#windows[@]}
echo "Windows: $nw"

################ write centers ################
> centers.dat
for c in "${windows[@]}"; do
  echo $c >> centers.dat
done

################################################
# LOOP OVER INCREASING END FRAME
################################################

start=2

for end in $(seq 51 50 501); do
  tag="2_${end}"
  echo "Processing NR=${start}..${end}"

  rm -f potall_${tag}.txt densityall_${tag}.txt densityall_${tag}.err
  rm -f *_us.dat

  ################ new file writing ###############
  for c in "${windows[@]}"; do
    awk -v s=$start -v e=$end 'NR>=s && NR<=e' \
        ../US-$c/$c"_us.dat" > $c"_us.dat"
  done

  ################ build potall.txt ################
  img=0
  for c in "${windows[@]}"; do

    awk -v img=$img -v n=$nw -v k1=$k_high -v k2=$k_low '
    NR==FNR {
        c[NR-1] = $1
        if ($1 <= 6.4)
            k[NR-1] = k1
        else
            k[NR-1] = k2
        next
    }
    NR > 1 {
        r = $8
        t = $1 * 2

        printf "%d %d %d", img, img, t
        for (i = 0; i < n; i++) {
            u = 0.5 * k[i] * (r - c[i]) * (r - c[i])
            printf " %.8f", u
        }
        printf "\n"
    }' centers.dat $c"_us.dat" >> potall_${tag}.txt

    img=$((img+1))
  done

  l=$(wc -l < potall_${tag}.txt)

  awk '{for(i=1;i<=NF;i++) if(i!=2) printf "%s ",$i; printf"\n"}' \
      potall_${tag}.txt | \
      $NPWHAM -w $nw -l $l -i 10000 -t $T \
      > densityall_${tag}.txt \
      2> densityall_${tag}.err

done

echo "ALL DONE"

cat densityall_2_51.err | tail -36 | awk '{print $3,$4}' > 1.dat
cat densityall_2_101.err | tail -36 | awk '{print $3,$4}' > 2.dat
cat densityall_2_151.err | tail -36 | awk '{print $3,$4}' > 3.dat
cat densityall_2_201.err | tail -36 | awk '{print $3,$4}' > 4.dat
cat densityall_2_251.err | tail -36 | awk '{print $3,$4}' > 5.dat
cat densityall_2_301.err | tail -36 | awk '{print $3,$4}' > 6.dat
cat densityall_2_351.err | tail -36 | awk '{print $3,$4}' > 7.dat
cat densityall_2_401.err | tail -36 | awk '{print $3,$4}' > 8.dat
cat densityall_2_451.err | tail -36 | awk '{print $3,$4}' > 9.dat
cat densityall_2_501.err | tail -36 | awk '{print $3,$4}' > 10.dat


paste 1.dat 2.dat 3.dat 4.dat 5.dat 6.dat 7.dat 8.dat 9.dat 10.dat | awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20}' > all.dat

awk '{
  sum=0; sumsq=0; n=0
  for (i=2; i<=7; i++) {
    sum += $i
    sumsq += $i*$i
    n++
  }
  avg = sum/n
  std = sqrt(sumsq/n - avg*avg)
  printf "%s %.6f %.6f\n", $1, avg, std
}' all.dat > fe.dat

rm 1.dat 2.dat 3.dat 4.dat 5.dat 6.dat 7.dat 8.dat 9.dat 10.dat all.dat 

paste centers.dat fe.dat > fe1.dat 
awk '{print$1,$3,$4}' fe1.dat > fe.dat
rm fe1.dat
