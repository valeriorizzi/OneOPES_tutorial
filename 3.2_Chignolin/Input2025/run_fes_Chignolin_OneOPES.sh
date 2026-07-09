
cd ../OneOPES_a/0
rm fes_stride*
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR.* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000 --min 0.015 --max 0.85;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_b/0
rm fes_stride*
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR.* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000 --min 0.015 --max 0.85;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_c/0
rm fes_stride*
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR.* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000 --min 0.015 --max 0.85;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_d/0
rm fes_stride*
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR.* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000 --min 0.015 --max 0.85;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES_e/0
rm fes_stride*
./../../FES_from_Reweighting.py --skip 20000 --sigma 0.01 --colvar COLVAR.* --cv rmsd_ca --bin 120 --temp 340 --deltaFat 0.2 --out fes_stride.dat --stride 20000 --min 0.015 --max 0.85;
grep 'DeltaF' fes_stride* | awk '{print -$4}' > deltaFl.dat
head -n 49 deltaFl.dat > deltaF.dat
cd ../../OneOPES

paste ../OneOPES_a/0/deltaF.dat ../OneOPES_b/0/deltaF.dat ../OneOPES_c/0/deltaF.dat ../OneOPES_d/0/deltaF.dat ../OneOPES_e/0/deltaF.dat > deltaF_list.dat;

awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*20000+20000, sum}' deltaF_list.dat > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF-1; sum3 = sqrt(sum2); print sum3}' deltaF_list.dat > temp2;
paste temp1 temp2 > deltaFall.dat;
rm temp*
