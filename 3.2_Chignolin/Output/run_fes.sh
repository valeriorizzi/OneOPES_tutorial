rm delta* fes_*

for i in a b c d e; do
    python3 ../../FES_from_Reweighting_multiT_funnel.py --skiprows 4000 --sigma 0.01 --bias opes.bias --colvar COLVAR.$i --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --min 0.015 --max 0.85 --blocks 3 --outfile fes_blocks_$i.dat
    python3 ../../FES_from_Reweighting_multiT_funnel.py --skiprows 4000 --sigma 0.01 --bias opes.bias --colvar COLVAR.$i --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --min 0.015 --max 0.85 --stride 500 --outfile fes_stride_$i.dat
    wait
    grep 'DeltaF ' fes_stride_$i* | awk '{print $4}' > deltaF_$i.dat
done

#wait

#append all the deltaF in one file
paste deltaF_* > deleteme;

#calculate average and stdev of all replicas in time
#Trypsin
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*500+4000, sum}' deleteme > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF; print sum2}' deleteme > temp2;
paste temp1 temp2 > deltaF_stride.dat;
sed -i '1i# time(ps) deltaF_stride(kJ/mol) ' deltaF_stride.dat
rm temp* deleteme
