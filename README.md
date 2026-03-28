# Human_IYD_Halogen_Bonding
Human IYD shows no halogen bonding in crystal structures, yet kinetics reveal strong halogen-dependent binding. Heavier halotyrosines bind tighter and associate faster than lighter halotyrosines, while koff changes little, suggesting halogens mainly boost kon via transient interactions. Static structures miss these dynamic effects. Cite: S. Karmakar and S. Mishra, J. Chem. Inf. Model., 2026 (under review).

<br>
<br>

The previous studies on Human IYD (both computational and experimental) showed that halogen bonding is absent in the active site of the enzyme. But, kinetic trends from experimental studies clearly highlighted a hidden role of halogen bonding in the substrate recognition pathway of the enzyme. Hence, we considered three substrates: I-Tyr (strong halogen and hydrogen bonding possible), Cl-Tyr (weak halogen and strong hydrogen bonding possible), and I-Phenol (strong halogen and weak hydrogen bonding possible). 

<br>
<br>

A. The parameter files were created in the same process as discussed in the Human_IYD_Dynamics (https://github.com/compchemsk/Human_IYD_Dynamics.git) repository. The parameter files for the three systems are given here in the Parameters directory. The input files for running classical MD simulations are similar to the files present in the Human_IYD_Dynamics directory.   

<br>
<br>

B. We first run a 1 microsec classical MD simulation for each system. From the simulation, we pick up a minimum free energy structure from the free energy landscape spanned by first two principal components. The detailed methodology and corresponding commands as well as codes are provided in the Human_IYD_Dynamics repository. At this stage, we are interested in isotropic steered MD (iSMD) simulations. For, iSMD simulations, a pre-equilibration step was performed. The input files are given in the iSMD directory. To run the simulation, follow the regular commands or commands are provided in the Human_IYD_Dynamics directory. 

<br>

The above simulation will generate a restart file named as Pre-Equilibration.rst7 which will be used in the next stage for pulling simulations. We will use different force constants for performing pulling simulations ranging from 0.1 to 4 kcal.mol-1.A-2. At each force constant value, we will run 100 pulling simulations using the code Pull_jobs.sh using the following command:

<br>
<br>
chmod +x Pull_jobs.sh
<br>
./Pull_jobs.sh

<br>
<br>

All the files required by the above code are present in the iSMD directory (change job.sh file as per your own computing node configurations). 

<br>
<br>
<br>

C. After the iSMD simulations are done, it is now time to check which force constant value is the Goldilocks limit for the ligand unbinding. Here, we use 6 different analyses to validate an appropriate force constant. For each force constant, we should have different directories (e.g., k0.1 to k4).  
<br>
1. Work distribution:
<br>
For work distribution accross different force constant, run the following codes in the sequence from the directory where k0.1 to k4 directories are present:
<br>
chmod +x Distribution.sh
<br>
./Distribution.sh
<br>
python work-dist-plot.py

<br>
A nice plot in python will appear showing work distribution. 

<br>
<br>

2. Free energy distribution:
<br>
To obtain free energy distribution, run the following codes in the same directory as before.
<br>
python free_energy_smd.py


<br>
<br>
<br>

3. Force distribution:
<br>
To obtain force distribution, run the following codes in the same directory as before.
<br>

for np in {0.1,0.25,0.5,1,2,3,4}
<br>
do
<br>
cp k_F.sh k${np}
<br>
cd k${np}
<br>
chmod +x k_F.sh
<br>
./k_F.sh
<br>
cd ../
<br>
done

<br>
<br>

After the above code, run:
<br>
python k_F_plot.py

<br>


<br>
<br>
<br>

4. PMF distribution:
<br>
To obtain PMF distribution, run the following codes in the same directory as before.
<br>

for np in {0.1,0.25,0.5,1,2,3,4}
<br>
do
<br>
cp k_PMF.sh k${np}
<br>
cd k${np}
<br>
chmod +x k_PMF.sh
<br>
./k_PMF.sh
<br>
cd ../
<br>
done

<br>
<br>

After the above code, run:
<br>
python k_PMF_plot.py

<br>

5. Motion of internal degrees of freedom:
<br>
To obtain free energy landscape spanned by internal degrees of freedom, run the following codes in the same directory as before. 
<br>

for np in {0.1,0.25,0.5,1,2,3,4}
<br>
do
<br>
cp int-coord.sh k${np}
<br>
cd k${np}
<br>
chmod +x int-coord.sh
<br>
sed -i "s/kxxx/k${np}/g" int-coord.sh
<br>
./int-coord.sh
<br>
cd ../
<br>
done

<br>
<br>

After the above code, go inside each directory (k0.1,k0.25,etc.) manually, and run:
<br>
chmod +x int-coord-fes.sh
<br>
./int-coord-fes.sh

<br>
