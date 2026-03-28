rm -rf Halogen_Bond_Evolution.dat Halogen_Bond_Evolution1.dat
vals=(3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8)

for ((i=0; i<${#vals[@]}-1; i++)) 

do

np=${vals[$i]}
next_np=${vals[$i+1]}

awk "/For trajectory ${np}/,/For trajectory ${next_np}/" halo.out | grep -v 'For trajectory' | grep -v 'Frame' | awk '{gsub(/^ +/,""); match($0,/Distance: ([0-9.]+)/,d); match($0,/Angle: ([0-9.]+)°/,a); print $1, d[1], a[1]}' | sort -k1,1 > ${np}.dat

awk -v np="$np" '{
    res=$1
    dist=$2
    ang=$3
    score=(ang/180)*(0.215/dist)

    sum[res]+=score
    sumsq[res]+=score*score
    count[res]++
}
END{
    printf "%-10s %-8s %-20s %-20s %-15s\n","Resname","Window","Mean_Geom_Quality_Score","Stddev_Geom_Quality_Score","Percent_Occupancy"
    for(r in count){
        occupancy=(count[r]/500)*100
        if(occupancy > 5){
            mean=sum[r]/count[r]
            std=sqrt((sumsq[r]/count[r]) - (mean*mean))
            printf "%-10s %-8s %-20.6f %-20.6f %-15.2f\n", r, np, mean, std, occupancy
        }
    }
}' ${np}.dat >> Halogen_Bond_Evolution1.dat
rm ${np}.dat
done
grep -v 'Resname' Halogen_Bond_Evolution1.dat > Halogen_Bond_Evolution.dat
rm Halogen_Bond_Evolution1.dat

awk '
NR==FNR {
    key = $1 FS $2
    $1=""; $2=""              # remove first two columns
    sub(/^  */,"")            # clean leading spaces
    val[key] = $0             # store remaining columns
    next
}
{
    key = $1 FS $2
    if (key in val)
        print $0, val[key]
}
' role-halo.dat Halogen_Bond_Evolution.dat > Halogen_Bond_Evolution2.dat
awk '{print$1,$2,$11,$10,$9/10,$3,$4,$5}' Halogen_Bond_Evolution2.dat > Halogen_Bond_Evolution3.dat
rm Halogen_Bond_Evolution.dat Halogen_Bond_Evolution2.dat

awk '{
    name=$1

    if(name=="FRA440"){
        $1="flavin"
    }
    else if(match(name, /^([A-Z]+)([0-9]+)$/, arr)){
        res=arr[1]
        num=arr[2]

        if(num < 219)
            $1 = res (num + 71)
        else if(num > 219)
            $1 = res (num - 149) "*"
    }

    print
}' Halogen_Bond_Evolution3.dat > Halogen_Bond_Evolution.dat

rm Halogen_Bond_Evolution3.dat
