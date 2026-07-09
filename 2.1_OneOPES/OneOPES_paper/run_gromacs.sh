#LIBTORCH=~/programs/plumed2-2.9/libtorch
#source ${LIBTORCH}/sourceme.sh;source /home/rizziv@farma.unige.ch/programs/plumed2-2.9/sourceme.sh;source /home/rizziv@farma.unige.ch/programs/gromacs-2022.3/install_single_cpu/bin/GMXRC
#export OMP_NUM_THREADS=1

create_tpr() {
  cd $i
  gmx_mpi grompp -f ../md.mdp -c ../config.gro -p ../topol.top -o prd.tpr
}

array=(`seq 0 7`)

rm -f */*.xtc */*.log */*.cpt */*.edr */COLVAR* */compres* */K* */delta* */fes* */bck* */#* */DeltaF*

for i in ${array[@]}
do
  create_tpr &
  #echo $i
done
wait

#10 ns simulation
mpirun -n 8 gmx_mpi mdrun -s prd.tpr -plumed plumed.dat -pme cpu -nb cpu -bonded cpu -pin on -pinoffset 16 -nsteps 5000000 -multidir 0 1 2 3 4 5 6 7 -hrex -replex 500 &
wait

cd 0
python3 ../../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR.* --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --stride 500 --out fes_phi.dat;
grep ' DeltaF ' fes_phi_* | awk '{print $4}' > deltaFstride.dat

