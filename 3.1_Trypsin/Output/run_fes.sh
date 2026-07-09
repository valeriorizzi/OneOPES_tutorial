rm delta* fes_*

for i in a b c d e; do
    python3 ../../FES_from_Reweighting_multiT_funnel.py --skiprows 2500 --sigma 0.01 --bias opes.bias --colvar COLVAR.$i --cv d1.z --bin 300 --temp 300 --min 0.2 --max 2.0 --rfunnel 0.2 --uat 1.8 --bat 0.8 --blocks 3 --outfile fes_blocks_$i.dat
    python3 ../../FES_from_Reweighting_multiT_funnel.py --skiprows 2500 --sigma 0.01 --bias opes.bias --colvar COLVAR.$i --cv d1.z --bin 300 --temp 300 --min 0.2 --max 2.0 --rfunnel 0.2 --uat 1.8 --bat 0.8 --stride 500 --outfile fes_stride_$i.dat
    wait
    grep 'DeltaF ' fes_stride_$i* | awk '{print $4}' > deltaF_$i.dat
done

#wait

#append all the deltaF in one file
paste deltaF_* > deleteme;

#calculate average and stdev of all replicas in time
#Trypsin
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*500+2500, sum}' deleteme > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF; print sum2}' deleteme > temp2;
paste temp1 temp2 > deltaF_stride.dat;
sed -i '1i# time(ps) deltaF_stride(kJ/mol) ' deltaF_stride.dat
rm temp* deleteme
