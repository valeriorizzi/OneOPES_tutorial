cd 0
gmx_mpi grompp -f ../production.mdp  -c ../B.pdb -p ../nosolv_GMX.top -o prd.tpr -r ../B.pdb -maxwarn 3
#repeat for all replicas

mpirun -n 8 gmx_mpi mdrun -deffnm prd -multidir 0 1 2 3 4 5 6 7 -replex 1000 -hrex -nsteps 125000000  -plumed plumed.dat -notunepme -pin on
