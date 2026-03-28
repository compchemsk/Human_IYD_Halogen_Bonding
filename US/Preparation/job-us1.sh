or np in {3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4}
do
mkdir US-${np}
cp hIYD-I-Tyr.parm7 US-${np}
cp rst7-input2/${np}-pull66.rst7 US-${np}
cp xx_us.in ${np}_us.in
sed -i "s/xx/${np}/g" ${np}_us.in
cp xx_us_dist1.RST ${np}_us_dist.RST
sed -i "s/xx/${np}/g" ${np}_us_dist.RST
cp ${np}_us.in ${np}_us_dist.RST US-${np}
rm ${np}_us.in ${np}_us_dist.RST
cd US-${np}
echo "pmemd.cuda -O -i ${np}_us.in -o ${np}_us.out -p hIYD-I-Tyr.parm7 -c ${np}-pull66.rst7 -r ${np}-us.rst7 -x ${np}-us.nc -inf ${np}-us.mdinfo" > ${np}-job-us.sh
chmod +x ${np}-job-us.sh
jsub ${np}-job-us.sh 1 gpu
cd ../
sleep 1s
done
