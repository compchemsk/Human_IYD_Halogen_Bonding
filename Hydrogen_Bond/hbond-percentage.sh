for np in 108 112 126 144 168 169 172 173 174 178 188 280 34 346 349 356 357 358 359 360 361 362 363 441 87 90 91 93 94 95 96 97 98 99
do
    awk -v val="$np" '$3 == val' hb-c1-ityrx.dat > ${np}-hb.dat

    np_values="3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8"
    output_file="${np}-hb-full.dat"

    # Print header
    echo -n "Resid" > "$output_file"
    for val in $np_values; do
        echo -n " $val" >> "$output_file"
    done
    echo >> "$output_file"

    # Generate matrix with AWK
    awk -v np_list="$np_values" '
    BEGIN {
        split(np_list, nps)
    }
    {
        val[$3][$1] = $4
        residues[$3] = 1
    }
    END {
        for (res in residues) {
            printf "%s", res
            for (i = 1; i <= length(nps); i++) {
                v = val[res][nps[i]]
                printf " %.2f", (v == "" ? 0 : v)
            }
            printf "\n"
        }
    }' ${np}-hb.dat >> "$output_file"

    rm ${np}-hb.dat
done

# Combine all into one file
cat *-hb-full.dat > hb-c1-it.dat
rm *-hb-full.dat
grep -v 'Resid' hb-c1-it.dat > hb-c1-itx.dat
rm hb-c1-it.dat 
mv hb-c1-itx.dat hb-c1-it.dat
