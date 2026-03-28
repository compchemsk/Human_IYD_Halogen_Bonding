for np in ARG93 ARG97 ASN171 ASN362 ASN89 FRA440 GLU363 GLY173 GLY95 HIS348 HIS359 HIS96 LEU188 LYS356 LYS357 LYS92 MET94 THR168 TRP98 TYR360 TYR361 VAL347 VAL358
do
    # Extract m and n from ${np}.dat file
    m=$(awk '{print $1}' ${np}.dat)
    n=$(awk '{print $2}' ${np}.dat)
    a=$((${n} + 1 ))
    # Copy necessary files to the target directory
    cp ../US-${m}/hIYD-I-Tyr.parm7 ../qm-halo/
    cp ../US-${m}/${m}-us.nc ../qm-halo/

    # Create the cpptraj input file
    echo "
parm hIYD-I-Tyr.parm7
trajin ${m}-us.nc ${a} ${a} 1
trajout ${a}.rst7 rst7
run
quit " > cp.in

    # Run cpptraj
    cpptraj cp.in

    # Calculate p value from np (extract number and increment)
    p=$(( ${np##*[!0-9]} + 1 ))

    # Extract the prefix (e.g., ARG, HIS, LYS) from np
    prefix=${np%%[0-9]*}

    # Set charge based on prefix
    case $prefix in
        ARG|LYS)     q=1 ;;
        GLU|ASP)     q=-1 ;;
        *)           q=0 ;;
    esac

    # Adjust q value
    q=$(( q - 1 ))

    # Create the minimization input file
    echo "
Initial min of our structure QMMM
 &cntrl
  imin=1, maxcyc=1,
  cut=8.0, ntb=1, ntc=2, ntf=2,
  ifqnt=1
 /
 &qmmm
  qmmask=':${p},442',
  qmcharge=${q},
  qm_theory='PM6', writepdb=1,
  qmshake=1,
  qm_ewald=1, qm_pme=1
 / " > min.in

    # Create the job script for sander
    echo "sander -O -i min.in -o min.out -p hIYD-I-Tyr.parm7 -c ${a}.rst7 -r min.rst7" > job.sh
    chmod +x job.sh

    # Run the job
    ./job.sh

    # Clean up temporary files
    rm min.in min.out hIYD-I-Tyr.parm7 ${a}.rst7 cp.in ${m}-us.nc mdinfo

    # Rename the PDB file and remove the first line
    mv qmmm_region.pdb ${np}-qm.pdb
    sed -i '1d' ${np}-qm.pdb
    ed -s ${np}-qm.pdb <<< $'$-3d\nw'
   
    # Extract atom coordinates (columns 3, 6, 7, 8) from PDB and save to qm.dat
    awk '{print $3, $6, $7, $8}' ${np}-qm.pdb > qm.dat

    # Create the Gaussian input file
    echo "%nprocshared=40
%mem=80GB
%chk=${np}-nbo.chk

#p b3lyp/genecp pop=(savenbo) scf=xqc empiricaldispersion=gd3 scrf=(cpcm,solvent=water,read) 

Title Line

${q} 1 " > 1.dat

    # Add the light atoms and basis set (6-311++G(d,p))
    echo "C H N O 0
6-311++G(d,p)
****
I     0
S    2   1.00
      0.7242000             -2.9731048
      0.4653000              3.4827643
S    1   1.00
      0.1336000              1.0000000
P    2   1.00
      1.2900000             -0.2092377
      0.3180000              1.1035347
P    1   1.00
      0.1053000              1.0000000
P    1   1.00
      0.0308                 1.0
D    1   1.00
      0.2940                 1.0
****

I     0
I-ECP     3     46
f potential
  5
0      1.0715702             -0.0747621
1     44.1936028            -30.0811224
2     12.9367609            -75.3722721
2      3.1956412            -22.0563758
2      0.8589806             -1.6979585
s-f potential
  5
0    127.9202670              2.9380036
1     78.6211465             41.2471267
2     36.5146237            287.8680095
2      9.9065681            114.3758506
2      1.9420086             37.6547714
p-f potential
  5
0     13.0035304              2.2222630
1     76.0331404             39.4090831
2     24.1961684            177.4075002
2      6.4053433             77.9889462
2      1.5851786             25.7547641
d-f potential
  5
0     40.4278108              7.0524360
1     28.9084375             33.3041635
2     15.6268936            186.9453875
2      4.1442856             71.9688361
2      0.9377235              9.3630657

Eps=4.0


" > 2.dat

    # Combine all the files to create the final Gaussian input file
    cat 1.dat qm.dat 2.dat > ${np}-nbo.com

    # Create the job script for Gaussian calculation
    echo "g16 ${np}-nbo.com" >> c1-g16.sh

    # Clean up temporary files
    rm 1.dat qm.dat 2.dat ${np}-qm.pdb job.sh 
sleep 1s
done
sed -i 's/C H N O 0/C H N S O 0/g' MET*com
sed -i 's/C H N O 0/C H N S O 0/g' CYS*com
sed -i 's/C H N O 0/C H N P O 0/g' FRA*com
rm *dat
