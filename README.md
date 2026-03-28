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

