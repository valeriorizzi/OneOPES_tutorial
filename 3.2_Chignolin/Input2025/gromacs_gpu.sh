#!/bin/sh
#SBATCH --partition=private-gervasio-gpu
#SBATCH --time 96:00:00                  # maximum run time.
#SBATCH --gpus=1
#SBATCH --constraint="V7|V8"
#SBATCH --job-name chigcharmm             # this is a parameter to help you sort your job when listing it
#SBATCH --error jobname-error.e%j     # optional. By default a file slurm-{jobid}.out will be created
#SBATCH --output jobname-out.o%j      # optional. By default the error and output files are merged
#SBATCH --ntasks 8                    # equivalent to -n for mpi
#SBATCH --cpus-per-task 1             # openmp threads 
#SBATCH --nodes 1             
export OMP_NUM_THREADS=1

alias splumedNov24='module load GCC/12.3.0 && module load OpenMPI && module load CUDA/12.3.0 && module load Python && module load SciPy-bundle && source /srv/beegfs/scratch/shares/flg/programs/plumed2-2.9.1_fixghost_1124/sourceme.sh; source /srv/beegfs/scratch/shares/flg/programs/gromacs-2023/install_mpi_1124/bin/GMXRC && module load OpenMPI'

splumedNov24
srun gmx_mpi mdrun -deffnm prd -multidir 0 1 2 3 4 5 6 7 -replex 10000 -hrex -plumed plumed.dat -nsteps 125000000 -cpi 

