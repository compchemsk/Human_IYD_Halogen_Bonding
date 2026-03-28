for np in TYR360 ASN171 ASN89 GLU363 MET94 ARG93 LEU188 THR168 HIS348
do
formchk ${np}-nbo.chk > ${np}-nbo.fchk
echo " 
import os
import subprocess

# Define filenames
formchk_file = 'xxx-nbo.fchk'
multiwfn_input = '''\
20
1
3
3
0
0
2
2
3
4
5
8
-4
6
0
-5
6
0
7
-1
-10
100
2
1
mol.pdb
'''

# Write input file for Multiwfn
with open('multiwfn_input.txt', 'w') as f:
    f.write(multiwfn_input)

# Run Multiwfn with formchk and automated input
print('Running Multiwfn...')
subprocess.run(f'Multiwfn {formchk_file} < multiwfn_input.txt', shell=True)

# Check if output files are created
outputs = ['mol.pdb', 'CPs.pdb', 'func1.cub', 'func2.cub']
missing = [f for f in outputs if not os.path.exists(f)]
if missing:
    print('Warning: Missing files:', ', '.join(missing))
else:
    print('Multiwfn output files generated successfully.')

# Launch VMD with both visualization scripts
print('Launching VMD...')
subprocess.run('vmd_new -e combine.vmd', shell=True) " > AIM-NCI.py
sed -i "s/xxx/${np}/g" AIM-NCI.py
python AIM-NCI.py
rm func1.cub func2.cub CPs.pdb paths.pdb mol.pdb AIM-NCI.py CPprop.txt *fchk multiwfn_input.txt
done
