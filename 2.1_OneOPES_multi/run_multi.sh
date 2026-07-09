#!/bin/bash
export OMP_NUM_THREADS=1

rm -f */*.xtc */*.log */*.cpt */*.edr */COLVAR* */compres* */K* */delta* */fes* */bck* */#* */DeltaF* delta* */mdout.mdp

total_iterations=5

array=($(seq 0 7))

create_tpr() {
  cd "$1" || exit
  gmx_mpi grompp -f ../md.mdp -c ../config.gro -p ../topol.top -o prd.tpr
}

for X in $(seq 1 $total_iterations); do
  echo "--- Starting iteration ${X} ---"

  for i in "${array[@]}"; do
    create_tpr "$i" &
  done
  wait

  mpirun -n 8 --oversubscribe gmx_mpi mdrun -s prd.tpr -plumed plumed.dat -pme cpu -nb cpu -bonded cpu -pin on -pinoffset 16 -nsteps 5000000 -multidir 0 1 2 3 4 5 6 7 -hrex -replex 200 &
  sleep 2

  echo ""
  echo "================================================================"
  echo " >>> Currently running iteration ${X} out of ${total_iterations}... <<<"
  echo "================================================================"
  echo ""
  wait

  cd 0 || exit
  python3 ../../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR.0 --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --stride 500 --out "fes_phi_iter${X}.dat"
  python3 ../../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR.0 --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --out "fes_blocks_phi_iter${X}.dat" --blocks 3;

  grep ' DeltaF ' fes_phi_iter${X}* | awk '{print $4}' > "deltaFstride_iter${X}.dat"

  cd ..

  if [ -f "0/COLVAR.0" ]; then
    mv 0/COLVAR.0 "0/COLVAR.0_iter${X}"
  fi

done

echo "=========================================="
echo " All ${total_iterations} iterations successfully completed."
echo "=========================================="

#append all the deltaF in one file
paste 0/deltaFstride_iter* > deleteme;

#calculate average and stdev of all replicas in time
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*500+1000, sum}' deleteme > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF; print sum2}' deleteme > temp2;
paste temp1 temp2 > deltaF_stride.dat;
sed -i '1i# time(ps) deltaF_stride_average(kJ/mol) deltaF_stride_stdv(kJ/mol)' deltaF_stride.dat
rm temp* deleteme

#block deltaF
grep ' DeltaF ' */fes_blocks_phi_iter* | awk '{v[NR]=$4; sum+=$4} END {mean=sum/NR; for(i=1;i<=NR;i++) sum2+=(v[i]-mean)^2; print mean, sqrt(sum2/(NR-1))}' > deltaF_block.dat
sed -i '1i# deltaF_block_average(kJ/mol) deltaF_block_stdv(kJ/mol)' deltaF_block.dat

echo "Script finished!"





