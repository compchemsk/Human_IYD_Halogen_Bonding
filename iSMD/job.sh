#!/bin/bash
#SBATCH -J xxx.sh                   #Job name(--job-name)
#SBATCH -o %j-out-myjob          #Name of stdout output file(--output)
#SBATCH -e %j-err-myjob          #Name of stderr error file(--error)
#SBATCH -p gpu        #Queue (--partition) name 
#SBATCH --nodes=1                   #number of nodes should be 1 for serial)
#SBATCH -n 1                #(--cpus-per-task) Number of Threads
#SBATCH --gres=gpu:1               # request gpu card: it should be either 1 or 2
#SBATCH -t 24:00:00         # specifies walltime(--time maximum duration)of run
##SBATCH --mem=23000        # Memory per node specification is in MB. It is optional. 
##SBATCH --mail-user=mishra@chem.iitkgp.ac.in        # user email ID where job status info will be sent
##SBATCH --mail-type=ALL        # Send Mail for all type of event regarding the job
# Thread Reservation Section
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
             omp_threads=$SLURM_CPUS_PER_TASK
           else
             omp_threads=1
fi
# Load Application/Compiler/Library modules, for example
# module load apps/lammps/12.12.2018/intel
 module load apps/gromacs/2020/cpu 
module load compiler/cuda/10.2 compiler/openmpi/4.0.2 compiler/openmpi/3.1.0/gnu
export AMBERHOME=/home/sabyashachicy/amber18
export PATH=$PATH:$AMBERHOME/bin
export PATH=$PATH:/home/sabyashachicy/software/bison/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sabyashachicy/software/bison/lib
export CUDA_HOME=/home/opt_ohpc_pub/apps/cuda/cuda-10.2
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib:$CUDA_HOME/lib64
# name of the executable
g16root="/home/sabyashachicy"
GAUSS_SCRDIR="/scratch/19cy40038/gaussian"
export g16root GAUSS_SCRDIR
. $g16root/g16/bsd/g16.profile
exe=./xxx.sh   
export OMP_NUM_THREADS=$omp_threads
# Executables along with right Parameters and right input files
#mpirun -bootstrap slurm -n $SLURM_NTASKS $exe  <   mnpc1_O2_HS.com  > mnpc1_O2_HS.log 
#export PATH=/home/apps/openmpi-4.0.2/bin:$PATH 
#export LD_LIBRARY_PATH=/home/apps/openmpi-4.0.2/lib:$LD_LIBRARY_PATH
#export orcadir=/home/sabyashachicy/software/orca_4_2_1_linux_x86-64_shared_openmpi314
#export PATH=$orcadir:$PATH
#export LD_LIBRARY_PATH=$orcalib:$LD_LIBRARY_PATH
nohup $exe
