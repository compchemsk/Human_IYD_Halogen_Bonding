for np in {1..100}
do
mkdir Pull-${np}
cp pull.in job.sh COM_pull.RST Pre-Equilibration.rst7 hIYD-I-Tyr.parm7 Pull-${np}
cd Pull-${np}
sed -i "s/Pull_dist.dat/Pull_dist_${np}.dat/g" pull.in
echo "pmemd.cuda -O -i pull.in -o pull-${np}.out -p hIYD-I-Tyr.parm7 -c Pre-Equilibration.rst7 -ref Pre-Equilibration.rst7 -r pull-${np}.rst7 -x pull-${np}.nc -inf pull-${np}.mdinfo" > I-Tyr-pull-${np}.sh
chmod +x I-Tyr-pull-${np}.sh
sed -i "s/xxx/I-Tyr-pull-${np}/g" job.sh
sbatch job.sh
sleep 1s
cd ..
done
