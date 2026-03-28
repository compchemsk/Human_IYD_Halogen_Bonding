for np in {6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.8,10.8,11.8,12.8,13.8,14.8,15.8,16.8,17.8,18.8}
do
mkdir US-${np}
cp hIYD-I-Tyr.parm7 US-${np}
cp rst7-input2/${np}-pull66.rst7 US-${np}
cp xx_us.in ${np}_us.in
sed -i "s/xx/${np}/g" ${np}_us.in
cp xx_us_dist2.RST ${np}_us_dist.RST
sed -i "s/xx/${np}/g" ${np}_us_dist.RST
cp ${np}_us.in ${np}_us_dist.RST US-${np}
rm ${np}_us.in ${np}_us_dist.RST
cd US-${np}
echo "pmemd.cuda -O -i ${np}_us.in -o ${np}_us.out -p hIYD-I-Tyr.parm7 -c ${np}-pull66.rst7 -r ${np}-us.rst7 -x ${np}-us.nc -inf ${np}-us.mdinfo" > ${np}-job-us.sh
chmod +x ${np}-job-us.sh
jsub ${np}-job-us.sh 1 gpu #Modify this line as per your system configuration
cd ../
sleep 1s
done
