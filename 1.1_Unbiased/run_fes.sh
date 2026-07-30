rm fes_phi* fes_b*

python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --blocks 5 --out fes_blocks_phi.dat;
grep 'DeltaF ' fes_blocks_phi.dat
grep 'DeltaF ' fes_blocks_phi.dat > deltaF_blocks.dat

python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --stride 5000 --out "fes_phi.dat"
grep ' DeltaF ' fes_phi_* | awk '{print $4}' > deltaFtemp.dat
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*5000, sum}' deltaFtemp.dat > deltaF_stride.dat
sed -i '1i# time(ns) deltaF_stride(kJ/mol) ' deltaF_stride.dat
rm deltaFtemp.dat
