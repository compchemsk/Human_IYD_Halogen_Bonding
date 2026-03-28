for np in {0.1, 0.25, 0.5, 1, 2, 3, 4} #change as per your force constant values
do
cp work-dist-smd.sh k${np}
cd k${np}
sed -i "s/xxx/${np}/g" work-dist-smd.sh
chmod +x work-dist-smd.sh
./work-dist-smd.sh
cp work-dist-k${np}.dat ../
cd ../
done
paste work-dist-k0.1.dat work-dist-k0.25.dat work-dist-k0.5.dat work-dist-k1.dat work-dist-k2.dat work-dist-k3.dat work-dist-k4.dat > work-dist-smd.dat
