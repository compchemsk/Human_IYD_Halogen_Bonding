cat ../hb-c1-ityrx.dat | sort -k3,3 | awk '$4>5' > Hydrogen_Bond_Evolution1.dat #hb-c1-ityrx.dat file contains progress_coord resid resval percent_occupan.

awk '{
    num = $3

    if(num == 441){
        $3 = "flavin"
    }
    else if(num < 220){
        $3 = num + 70
    }
    else if(num > 220){
        $3 = (num - 150) "*"
    }

    print
}' Hydrogen_Bond_Evolution1.dat > Hydrogen_Bond_Evolution2.dat

awk '{print $1,$3,$4}' Hydrogen_Bond_Evolution2.dat > Hydrogen_Bond_Evolution.dat
rm Hydrogen_Bond_Evolution1.dat Hydrogen_Bond_Evolution2.dat
