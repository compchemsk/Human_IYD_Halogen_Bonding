#target.sh will provide the values which will be included in this rst7-input.sh, inside the for loop.
rm -rf rst7-input
mkdir rst7-input
for np in {34,378,621,599,972,1013,1064,1070,1072,1088,1129,1159,1171,1182,1184,1189,1187,1185,1190,1188,1191,1192,1193}
do
echo "
parm hIYD-I-Tyr.parm7
trajin pull-66.nc ${np} ${np} 1 #Randomly selected trajectory from the cluster.
trajout rst7-input/${np}-pull66.rst7 rst7
run
quit " > cpptraj-${np}.in
cpptraj cpptraj-${np}.in
rm cpptraj-${np}.in
done
