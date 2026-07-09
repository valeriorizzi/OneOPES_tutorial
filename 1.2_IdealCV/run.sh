#!/bin/bash
export OMP_NUM_THREADS=1

mpirun -n 1 gmx_mpi mdrun -deffnm prd -plumed plumed.dat -pin on -ntomp 1 -pme cpu -nb cpu -bonded cpu -nsteps 5000000 &
wait

rm fes_*

python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --blocks 5 --out fes_blocks_phi.dat;
grep 'DeltaF ' fes_blocks_phi.dat
grep 'DeltaF ' fes_blocks_phi.dat > deltaF_blocks.dat

python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --stride 500 --out "fes_phi.dat"
grep ' DeltaF ' fes_phi_* | awk '{print $4}' > deltaFtemp.dat
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*500+1000, sum}' deltaFtemp.dat > deltaF_stride.dat
sed -i '1i# time(ps) deltaF_stride(kJ/mol) ' deltaF_stride.dat
rm deltaFtemp.dat
