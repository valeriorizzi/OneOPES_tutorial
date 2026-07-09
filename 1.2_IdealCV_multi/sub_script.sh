#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o /ELS/els9/users/ngmms2026/work/$USER/sge/out
#$ -e /ELS/els9/users/ngmms2026/work/$USER/sge/err
####$ -N threaded

echo "Running on host:"
hostname
source /ELS/els9/users/ngmms2026/ngmm2026_env.sh
module purge
#module load gromacs-2023.5
mopdule load gromacs-2024.2_plumed
module list 2>&1
echo "Allocated slots: $NSLOTS"
#export OMP_NUM_THREADS=$NSLOTS
export OMP_NUM_THREADS=1
#gmx_mpi -version

rm -f */*.xtc */*.log */*.cpt */*.edr */COLVAR* */compres* */K* */delta* */fes* */bck* */#* */DeltaF*
mpirun -n 4 gmx_mpi mdrun -deffnm prd -plumed ../plumed.dat -pin on -pinoffset 16 -ntomp 1 -pme cpu -nb cpu -bonded cpu -nsteps 5000000 -multidir 0 1 2 3

#qrsh -l seixi=1 -pe mpi 4 sub_script.sh
#qsub -l seixi=1 -pe seixi.32t mpi 4 sub_script.sh
#qsub -l seixi=1 -pe openmpi.4t 4 sub_script.sh
