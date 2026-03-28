#awk '{printf "\"%s\",", $0}' halo-amino.dat | sed 's/,$//'
#!/bin/bash

for np in ARG93 ARG97 ASN171 ASN362 ASN89 FRA440 GLU363 GLY173 GLY95 HIS348 HIS359 HIS96 LEU188 LYS356 LYS357 LYS92 MET94 THR168 TRP98 TYR360 TYR361 VAL347 VAL358
do
    grep "${np}" ../halo.out > 1.dat
    grep -oP 'Angle:\s*\K[0-9.]+' 1.dat > x.dat
    grep -oP 'Distance:\s*\K[0-9.]+' 1.dat > y.dat
    paste y.dat x.dat > ${np}.dat
    rm y.dat x.dat 1.dat

    # Sort by angle (column 2), take the max
    sort -nk2 ${np}.dat | tail -1 > ${np}-f.dat
    rm ${np}.dat

    # Read distance and angle
    read target_dist target_angle < ${np}-f.dat

    # Format
    dist_pattern=$(printf "%.3f" "$target_dist")
    angle_pattern=$(printf "%.2f" "$target_angle")

    trajectory=""
    frame=""
    found=0

    while IFS= read -r line; do
        if [[ $line =~ For\ trajectory\ ([0-9]+\.[0-9]+) ]]; then
            trajectory="${BASH_REMATCH[1]}"
        fi

        if [[ $line =~ Frame\ ([0-9]+): ]]; then
            frame="${BASH_REMATCH[1]}"
        fi

        if [[ $line =~ Distance:\ ${dist_pattern}\ nm ]] && [[ $line =~ Angle:\ ${angle_pattern}° ]]; then
            echo -e "${trajectory}\t${frame}\t${target_dist}\t${target_angle}" > ${np}-f.dat
            found=1
            break
        fi
    done < ../halo.out

    if [[ $found -eq 0 ]]; then
        echo "No matching frame found for ${np} with Distance=$dist_pattern and Angle=$angle_pattern"
    else
        mv ${np}-f.dat ${np}.dat
    fi
done
