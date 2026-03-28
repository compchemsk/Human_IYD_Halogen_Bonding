#!/bin/bash
#SBATCH -J jobname                    #Job name(--job-name)
#SBATCH -o %j-out-myjob          #Name of stdout output file(--output)
#SBATCH -e %j-err-myjob          #Name of stderr error file(--error)
#SBATCH -p medium        #Queue (--partition) name
#SBATCH --nodes=1                   #number of nodes should be 1 for serial)
#SBATCH -n 40                #(--cpus-per-task) Number of Threads
#SBATCH -t 72:00:00
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
export AMBERHOME=/home/sabyashachicy/amber18
export PATH=$PATH:$AMBERHOME/bin
export PATH=$PATH:/home/sabyashachicy/software/bison/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sabyashachicy/software/bison/lib
export CUDA_HOME=/home/opt_ohpc_pub/apps/cuda/cuda-10.2
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib:$CUDA_HOME/lib64
# name of the executable
g16root="/home/sabyashachicy"
GAUSS_SCRDIR="/scratch/sabyashachicy/gaussian"
export g16root GAUSS_SCRDIR
. $g16root/g16/bsd/g16.profile
export OMP_NUM_THREADS=$omp_threads
# Executables along with right Parameters and right input files
#mpirun -bootstrap slurm -n $SLURM_NTASKS $exe  <   mnpc1_O2_HS.com  > mnpc1_O2_HS.log
 
g16 FRA440-nbo.com
g16 ARG93-nbo.com
g16 ARG97-nbo.com
g16 ASN171-nbo.com
g16 ASN362-nbo.com
g16 ASN89-nbo.com
g16 FRA440-nbo.com
g16 GLU363-nbo.com
g16 GLY173-nbo.com
g16 GLY95-nbo.com
g16 HIS348-nbo.com
g16 HIS359-nbo.com
g16 HIS96-nbo.com
g16 LEU188-nbo.com
g16 LYS356-nbo.com
g16 LYS357-nbo.com
g16 LYS92-nbo.com
g16 MET94-nbo.com
g16 THR168-nbo.com
g16 TRP98-nbo.com
g16 TYR360-nbo.com
g16 TYR361-nbo.com
g16 VAL347-nbo.com
g16 VAL358-nbo.com
