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
<br>

6. During unbinding of a substrate from the active site, e.g., I-Tyr, crucial hydrogen bonding interactions are broken. To monitor the sequence of events, run the following code, after going inside each directory (k0.1,k0.25,etc.) manually:

<br>
<br>
chmod +x seq-of-events.sh
<br>
./seq-of-events.sh

<br>
<br>
Same analysis were repeated for Cl-Tyr and I-Phenol to obtain the optimal force constant values in each case. All these codes are available in the Goldilocks_Limit directory. 


<br>
<br>
<br>

D. After running non-equilibrium simulations (iSMD), we can extract kinetic information from them using Dudko-Hummer-Szabo (DHS) model. The codes are given in the DHS directory. To run the codes, the user needs 4 inputs as: rupture force values (from iSMD simulations), corresponding standard deviation values, substrate unbinding time (from simulation trajectories), and corresponding standard deviation values.   
<br>
<br>
<br>

E. Now, we have an optimal force constant value as 0.5 kcal/mol.A^2. We also have 100 unbinding trajectories using the same force constant. Hence, at this stage, we want to see how many different pathways emerge for substrate unbinding. To classify distinct substrate unbinding pathways explored in iSMD simulations, we employed Dynamic Time Warping (DTW) to measure similarity between trajectories that may differ in length and progression along the reaction coordinate. Each simulation was represented using two one-dimensional descriptors: substrate RMSD and an internal coordinate (center-of-mass distance between flavin and the substrate). DTW distances were computed separately for these descriptors and then averaged to obtain a composite similarity measure between trajectory pairs. The resulting pairwise distance matrix was subjected to agglomerative hierarchical clustering using Ward’s linkage criterion. The most populated cluster was identified as the dominant unbinding pathway, and a representative trajectory from this cluster was selected as the starting structure for subsequent umbrella sampling simulations. All the codes and inputs are available in the DTW directory. 

<br>
First run:
<br>
chmod +x int-dist-full.sh
<br>
./int-dist-full.sh

<br>
chmod +x rmsd_full.sh
<br>
./rmsd_full.sh

<br>
<br>
Finally run:
<br>
pip install fastdtw scipy seaborn matplotlib numpy
<br>
python DTW_Clustering.py

<br>
<br>
<br>
Same analysis were performed for Cl-Tyr and I-Phenol systems.


<br>
<br>
<br>

F. Umbrella sampling (US) simulations were performed across the progress coordinate using multiple force constants to ensure adequate sampling. Initial tests revealed that low force constants led to poor window overlap, while high values caused overbias. Based on this, an optimal scheme was adopted: a force constant of 10 kcal·mol⁻¹·Å⁻² with finer spacing (0.2 Å) at early stages, and 2.5 kcal·mol⁻¹·Å⁻² with wider spacing (1 Å) at larger distances (Optimum values can be identified by varrying the force constant, window size, and extent of simulation). Each window was then simulated for 25 ns, and PMFs were computed using non-parametric reweighting, reducing binning artefacts and yielding smoother free energy profiles. The input files and codes are given in the US directory.

<br>
<br>
Run the files in this sequence as below:
1. Create dis.dat file using cpptraj from parameter and trajectory file.
<br>
2. Run: a) target.sh, b) rst7-input.sh, c) name-change.sh, d) job-us-run.sh (to run use: chmod +x filename.sh and then ./filename.sh )  (All codes are given in the US/Preparation directory.)

<br>
<br>
<br>

Once, umbrella sampling simulations are done, we will use non-parametric reweighting scheme (npwham) (Please Cite: https://doi.org/10.1038/s43588-022-00389-9) to construct the PMF profiles. Run the npwham-run.sh code to get the PMF files.