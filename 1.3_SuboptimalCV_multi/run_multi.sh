#!/bin/bash
export OMP_NUM_THREADS=1

#gmx_mpi grompp -f ../md.mdp -c ../config.gro -p ../topol.top -o prd.tpr

# --- Configuration Variables ---
TIME_NS=10
DT=0.002

NSTEPS=$(awk -v t=$TIME_NS -v dt=$DT 'BEGIN { printf "%.0f", (t * 1000) / dt }')

folders=("0" "1" "2" "3")
file_list=""

rm -f */*.xtc */*.log */*.cpt */*.edr */COLVAR* */compres* */K* */delta* */fes* */bck* */#* */DeltaF* delta* */mdout.mdp

########## running calculation, to be commented if one wants to only analyse results
rm -f */*.xtc */*.log */*.cpt */*.edr */COLVAR* */compres* */K* */delta* */fes* */bck* */#* */DeltaF*
mpirun -n 4 --oversubscribe gmx_mpi mdrun -deffnm prd -plumed ../plumed.dat -pin on -pinoffset 16 -ntomp 1 -pme cpu -nb cpu -bonded cpu -nsteps $NSTEPS -multidir 0 1 2 3 &
wait
##########

# --- Per-Folder Analysis Loop ---
for dir in "${folders[@]}"; do
    echo "--- Processing folder: $dir ---"
    cd "$dir"

    (
        python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR.* --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --stride 500 --out fes_phi.dat;
        grep ' DeltaF ' fes_phi_* | awk '{print $4}' > deltaFstride.dat
        python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR.* --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --out fes_blocks_phi.dat --blocks 3;
    ) &

    file_path="${dir}/deltaFstride.dat"
    file_list+=" $file_path"
    cd ..
done

echo "Waiting for all background FES jobs to complete..."
wait
echo "All FES calculations are completed."

#stride deltaF
paste $file_list > deltaF_list2.dat;
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*500+1000, sum}' deltaF_list2.dat > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF-1; sum3 = sqrt(sum2); print sum3}' deltaF_list2.dat > temp2;
paste temp1 temp2 > deltaF_stride.dat;
sed -i '1i# time(ps) deltaF_stride_average(kJ/mol) deltaF_stride_stdv(kJ/mol)' deltaF_stride.dat
rm temp* deltaF_list2.dat

#block deltaF
grep ' DeltaF ' */fes_blocks_phi.dat | awk '{v[NR]=$4; sum+=$4} END {mean=sum/NR; for(i=1;i<=NR;i++) sum2+=(v[i]-mean)^2; print mean, sqrt(sum2/(NR-1))}' > deltaF_block.dat
sed -i '1i# deltaF_block_average(kJ/mol) deltaF_block_stdv(kJ/mol)' deltaF_block.dat

echo "Script finished!"
