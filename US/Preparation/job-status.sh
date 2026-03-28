for np in {3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.8,10.8,11.8,12.8,13.8,14.8,15.8,16.8,17.8,18.8}
do
cd US-${np}
rm job-stat.dat
cat *_us.out | tail -1 >> ../job-stat.dat
cd ../
done
cat job-stat.dat -n
rm job-stat.dat
