for i in {1..100}
do
cp 2d.py 2d_boltzman.f Pull-${i}
cd Pull-${i}
sed -i "s/xxx/internal_coord_k0.1_${i}.dat/g" 2D_Boltzman.f #change k0.1 to whatever directory, you are currently in
gfortran 2D_Boltzman.f
./a.out > kB-${i}.txt
sed -i "s/xxx/kB-${i}.txt/g" 2D.py
python 2D-fes.py
cd ../
done
