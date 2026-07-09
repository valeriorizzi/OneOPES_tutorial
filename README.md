# Tutorial: Mastering Enhanced Sampling with OneOPES

<div align="center">

<img width="350" height="350" alt="ezgif-49494ac320700ea2" src="https://github.com/obzehn/renders/blob/main/gifs/video_hostguest.gif?raw=true" />

*Figure 1: OneOPES logo and a host-guest system bound state dynamics. Image created and rendered by [Nicola Piasentin](https://github.com/obzehn/renders)*

</div>

This tutorial provides a technical guide to implementing **OneOPES** ([JCTC 2023](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00254)), a replica-exchange strategy developed by Valerio Rizzi, Simone Aureli, Narjes Ansari and Francesco Luigi Gervasio. OneOPES combines enhanced sampling techniques **OPES Explore** ([JCTC 2022](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00152)) and **OPES MultiThermal** ([PRX 2020](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034)), recently developed by Michele Invernizzi, Pablo Piaggi and Michele Parrinello. The input files from the OneOPES paper are all available on ([PLUMED NEST](https://www.plumed-nest.org/eggs/23/011/)) and ([Github](https://github.com/valeriorizzi/OneOPES/)).

OneOPES is designed to consistently deliver valid free energy results even when the accelerated Collective Variables (CVs) are suboptimal, significantly reducing the human effort usually required for fine-tuning reaction coordinates. The question of wether an enhanced sampling simulation has reached convergence is a crucial one. In fact, convergence is surprisingly hard to define unequivocally. When a reference is available, it is rather straightforward to assess convergence by direclty comparing the estimated free energy and with the reference. However, computational references are often only available in toy model systems that have an analytical solution or that can be simulated for a sufficiently long time. When dealing with open areas of research, such a reference is often unavailable. Establishing criteria to evaluate the quality of an enhanced sampling simulation is a paramount topic that we will discuss in the tutorial.

In OneOPES, we typically have a set of replicas that range in exploratory power. While main replica 0 contains only a single bias on primary CVs, higher replicas explicitly accelerate additional degrees of freedom. Therefore an accelerated exploration of the phase space is achieved in higher replicas. The visited configurations percolate down to replica 0 through replica exchange. Replica 0 thus gains in sampling quality and is used to calculate thermodynamic properties such as the free energy. A more detailed discussion on OneOPES is present in Section 2.

**Disclaimer**: OneOPES is not _magic_ and does not encourage the user to be sloppy and not spend time to improve their CVs whenever possible. On the contrary, like any other CV-based enhanced sampling method, the quality of its results improves when the CV quality improves. However, in large, complicated systems such as those found in biophysics, users often reach a point where any improvement to the CV comes at a high cost, i.e. it may be very difficult, time consuming and often system dependent. These complex systems are the true target of the OneOPES strategy. Some of its applications are listed in Section 3 at the end of the tutorial.

<img width="1302" height="550" alt="oneopes_scheme" src="https://github.com/user-attachments/assets/862345a0-4af0-40db-aefb-7b41b2116f9a" />

*Figure 2: Schematic representation of the OneOPES replica exchange method ([JCTC 2023](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00254))*

---

## Prerequisites

Ensure you have the following tools installed and configured:

*   **MD Engine:** GROMACS 2022.5 / 2023 patched with PLUMED 2.9.1
*   **Parallel Computing:** An MPI-enabled environment to run 8 replicas simultaneously.
*   **Analysis Tools:** Python 3, Gnuplot and GNU Awk for processing `COLVAR` files, performing reweighting and extracting thermodynamic properties.

The simulations are not computationally expensive. You can run them locally on your workstation or on an HPC partition, ideally with an interactive job.

Beginner users are encouraged to gain confidence with the syntax of PLUMED and the basics of statistics by following the tutorials [PLUMED Masterclass 21.1](https://www.plumed-tutorials.org/lessons/21/001/data/NAVIGATION.html) and [PLUMED Masterclass 21.2](https://www.plumed-tutorials.org/lessons/21/002/data/NAVIGATION.html).

---

For the NGMM2026 students, after logging into the HPC facility. Create an interactive job with
```bash
qrsh -l seixi=1 -pe mpi 2
```

and then in the terminal run
```bash
source /ELS/els9/users/ngmms2026/ngmm2026_env.sh
module purge
module load python-3.13.1
module load gromacs-2023-plumed-2.9.1-MPI
```

---

## Step 1: Understanding CV Quality with Alanine Dipeptide

Alanine dipeptide (Ala2) is an ideal benchmark system for a tutorial on enhanced sampling. This is particularly true in our case where we will discuss how CV quality dictates convergence. It is customary to define its conformation using the dihedral angles $\phi$ and $\psi$ phase space. The expected $\Delta F$ between its two main basins is $8.9 \pm 0.1$ kJ/mol ([JCTC 2023](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00254)).

<img width="1271" height="535" alt="ala2_valsson2016" src="https://github.com/user-attachments/assets/c5aac508-c8ba-4dba-a6ae-5e8405f96687" />

*Figure 3: Visual depiction of Alanine Dipeptide, its two metastable states and a reference 2-dimensional FES $`F(\phi, \psi)`$ in $`k_B T`$ ([AnnRev 2016](http://www.annualreviews.org/doi/10.1146/annurev-physchem-040215-112229)). Pay attention to the apparent kinetic barrier of about 20 $`k_B T`$ between the two metastable states. This barrier would make spontanoeus thermal transitions between them _rare events_, as we will see in the Unbiased trajectory presented in Section 1.1.*

**Another Disclaimer**: Ala2 is a useful playing ground to introduce important concepts, such as comparing different CVs, building a free energy surface, and calculating its error. However, it cannot be relied upon as a truly probing system. It is effectively a toy model that should not illude users to draw any solid conclusion about the general capabilities of an enahnced sampling method or of a CV building strategy on complex systems. Today, Ala2 is way too simple to represent a challenge for modern methods. Explicit tests on complex systems are needed to properly test a method. Such example exist for OneOPES and they are the true target of the strategy, but they are too expensive to run in a tutorial. We include in Section 3 a list of references to publications that use OneOPES and add the input and output of two examples: one of protein-ligand binding and one of miniprotein folding. In Section 2.2 we provide some parameter choice reccomentation and in the Troubleshooting section we point out possible solutions to common problems.

### 1.1 Unbiased reference

A 87 microseconds unbiased trajectory is available on [Zenodo](https://zenodo.org/records/7323535) as a reference.

1.1.1  **Visualise the reference** 

Enter folder 1.1_Unbiased where you will find a downsampled version of the reference COLVAR file. Open Gnuplot in the terminal and run
```bash
p 'COLVAR' u 1:2
```
to visualise the dynamics of $\phi$. 

    
In Gnuplot, type
```bash
p 'COLVAR' u 1:3:2 pal
```
to visualise the dynamics of $\psi$ coloured by $\phi$.


<blockquote>

[!NOTE]
    
<details>
<summary>❓🤔 Does $\psi$  distinguish well the Ala2 basins? How does it compare to $\phi$ </summary>
<br>
There is a strong overlap between the two metastable basins in $\psi$, while $\phi$ distinguishes them well. While there are a number of other factors that determine the quality of a CV to be accelerated in enhanced sampling, if a CV does not distinguish the free energy basins and the transition state region, it has no hope to be a good quality CV. See also the FES in Fig. 3.
</details>
    
</blockquote>

1.1.2  **Calculate the Free energy surface and the $\Delta F$** 

In the the main folder, you can find the [FES_from_Reweighting_multiT_funnel.py](https://github.com/obzehn/multithermal_fes/blob/main/FES_from_Reweighting_multiT_funnel.py) script. It is a modified version by Nicola Piasentin of the original [reweighting script](https://github.com/invemichele/opes/blob/master/postprocessing/FES_from_Reweighting.py) by Michele Invernizzi ([JPCL 2020](https://pubs.acs.org/doi/10.1021/acs.jpclett.0c00497)). Run the script
```bash
bash run_fes.sh
```
to calculate the 1D free energy surface (FES) over $\phi$ and the $\Delta F$ between the metastable basins. 

With this command
```bash
python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --blocks 5 --out fes_blocks_phi.dat;
grep 'DeltaF ' fes_blocks_phi.dat
```
we run a block average and extract the corresponding $\Delta F$ value.

And you can visualise the 1D FES with Gnuplot
```bash
p[][-2:100] 'fes_blocks_phi.dat' w e
``` 

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 How does the FES estimation changes if you vary the number of blocks? Try with the flag --blocks 50, for example. </summary>
<br>   
Each block must contain enough sampling of the two basins. If that is not the case, i.e. when you use a a too large number of blocks, the quality of the FES estimation dramatically worsens.
</details>    
</blockquote>

Another way to calculate the FES and the $\Delta F$ is to do so sequentially by using the option --stride
```bash
python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --stride 5000 --out "fes_phi.dat"
grep ' DeltaF ' fes_phi_* | awk '{print $4}' > deltaFtemp.dat
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*5000, sum}' deltaFtemp.dat > deltaF_stride.dat
```
This shows the progression of $\Delta F$ as a function of simulation time, but does not offer an error estimation.

### 1.2 The Ideal Case Scenario, biasing CV $\phi$

The angle **$\phi$** is a known high-quality CV for Ala2 as it captures well the two main basins and their transition state.

1.2.1  **Run a single OPES Explore simulation:** 

Enter folder 1.2_IdealCV and run command
```bash
bash run.sh
```  
You will perform a standard OPES Explore simulation biasing $\phi$ for 10 ns. Here we use a `BARRIER` of 30 kJ/mol and a `PACE` of 2000 steps. The script will run the simulation with this GROMACS command 
```bash
mpirun -n 1 gmx_mpi mdrun -deffnm prd -plumed plumed.dat -pin on -ntomp 1 -pme cpu -nb cpu -bonded cpu -nsteps 5000000 &
```  

The corresponding plumed.dat script is 

``` plumed
MOLINFO STRUCTURE=input.ala2.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
omega: TORSION ATOMS=@omega-2
theta: TORSION ATOMS=6,5,7,9
ene: ENERGY

OPES_METAD_EXPLORE ...
  LABEL=opes
  ARG=phi
  FILE=Kernels.data
  STATE_RFILE=compressed.Kernels
  STATE_WFILE=compressed.Kernels
  PACE=2000
  BARRIER=30
... OPES_METAD_EXPLORE

PRINT STRIDE=500 FILE=COLVAR ARG=opes.bias,phi,psi,omega,theta,ene FMT=%7.4f

ENDPLUMED
```

that is followed by a FES estimation, using the stride and the block approaches
```bash
python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --blocks 5 --out fes_blocks_phi.dat;
grep 'DeltaF ' fes_blocks_phi.dat > deltaF_blocks.dat

python3 ../FES_from_Reweighting_multiT_funnel.py --sigma 0.05 --colvar COLVAR --cv phi --bin 150 --temp 300 --min -3 --max 2 --deltaFat -0.2 --skip 1000 --stride 500 --out "fes_phi.dat"
grep ' DeltaF ' fes_phi_* | awk '{print $4}' > deltaFtemp.dat
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*500+1000, sum}' deltaFtemp.dat > deltaF_stride.dat
```

Here we are dealing with an enhanced sampling simulation that iteratively builds an external potential over the CVs to accelerate sampling. Note two key differences compared to the unbiased case. 
First, each frame has a statistical weight determined by the external potential contained in the column opes.bias. That is taken into account by the instruction --bias. This instruction is not compulsory in the script as, by default, every column *.bias existing in COLVAR is read, summed over and used for reweighting.
Second, the initial portion of the trajectory when the external bias is being built is out-of-equilibrium ([AnnRev 2016](http://www.annualreviews.org/doi/10.1146/annurev-physchem-040215-112229)) and must be discarded from the reweight procedure with keyword --skip. Here, we propose to exclude a standard amount, 10% of the trajectory. 

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 Do the block and the stride $\Delta F$ estimates agree with each other? Does it agree with the expected value of $8.9 \pm 0.1$ kJ/mol? </summary>
<br>   
Yes, the block and stride estimates do agree with each other and this is a reassuring sign that the simulation converged well. We also agree with the expected value that is another confirmation that the simulation converged well.
</details>    
</blockquote>   

<blockquote>    

[!NOTE]
<details>
<summary>❓🤔 Change the number of blocks from 5 to 50 as you did in the unbiased example. How does the $\Delta F$ estimation compare now? </summary>
<br>   
This simulation is much more robust to a change in the number of blocks, as the simulation has a high quality sampling of the two basins. However, the estimate can still worsen if the blocks' length is too short. A solution is either to prolong the simulation or to choose a more suitable block number.
</details>    
</blockquote>   

In Gnuplot, run
```bash
p 'COLVAR' u 1:3:2 pal        
``` 
to see the dynamics of $\phi$ colour by the opes.bias. 
    
> [!IMPORTANT]
> <span id="transitions-discussion">Many transitions are observed, which is a good sign, but is not sufficient for building a converged free energy. Counting the number of back and forth transition is in fact a deceiving criterium for convergence. It is true that back and forth transitions are necessary for building a free energy as one must sample well the relevant phase space. However, if those transitions are due to an out-of-equilibrium bias potential that strongly changes in time, the reweight procudere would be noisy, ineffective and overall ill-defined. One must have back and forth transitions _and_ a quasi-static underlying bias potential to build a converged FES. </span>

To check this, run
```bash
p 'COLVAR' u 1:2:3 pal        
``` 
to plot the dynamics of opes.bias coloured by the $\phi$ value. For an ideal reweighting, at long times one must have rather constant bias values for analogous CV values. In the plot, this is visible as constant colour bands at different levels. This is a reassuring sign that indicates the good quality of this simulation. 

If you focus on the very beginning of the simulation such as the first 1000 ps
```bash
p[:1000] 'COLVAR' u 1:2:3 pal        
``` 
you can see spikes at constant colour or CV value, indicating that the bias is out-of-equilibrium and that portion of trajectory is better to be skipped in the reweighting procedure. 

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 Compare the estimated $\Delta F$ with the one from the unbiased trajectory in Section 1.1. How long did it take to run the OPES Explore simulation? Which way do you think is more efficient? </summary>
<br>
The two estimates are comparable, but the enhanced sampling one is more accurate. In the unbiased trajectory, there are some back and forth transitions but not a large enough quantity to build a robust statistics. Given the long length of the simulation and its computational cost, this is far from an efficient method. To put the cost in perspective, it took about 5 days to run the unbiased simulations, while the OPES Explore simulation should have taken you less than a minute. Please note that this massive acceleration is due to Ala2 being small and the CV being optimal. Nowadays, in more realistic systems, unbiased simulations can compete in efficiency with enhanced sampling simulations, depending on a number of factors, including system size, computational resources, computational method and, especially, CV quality.
</details>

</blockquote>

  
1.2.2  **A more robust estimate of the error, run multiple independent simulations:**     

Reliably estimating errors from a single simulation poses a certain degree of danger, as the trajectory is biased by the same external potential so it is self-correlated to an extent. Independent trajectories build instead uncorrelated independent external potentials and they represent a stronger convergence test, thus providing a more robust error estimate. 

Change folder and enter 1.2_idealCV_multi. Here you will find 4 folders (0 1 2 3), each containing analogous starting tpr files that start from the same initial configuration with different initial velocities. By running
```bash
bash run_multi.sh
``` 
you will perform 4 independent simulations and calculate the FES and $\Delta F$ with blocks and stride. In the parent folder, file deltaF_stride.dat contains the mean and standard deviation of $\Delta F$ in time over the independent replicas. In each simulation, $\Delta F$ is calculated at a given time and we provide mean and standard deviation of those values. File deltaF_blocks.dat contains the same quantities calculated as a block average in each independent simulation. The script run_average_fes.sh will generate a 1D FES that calculates the average and stdev of each bin of each independent replicas.

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 Do the multiple simulation results agree with the previous single simulation ones? </summary>
<br>   
The two results should be in good agreement, which is a further indication of converged simulations. Performing independent simulation is a robust way to estimate errors and to check whether the simulations are truly converged. 
</details>   
</blockquote>   

### 1.3 The Bad Case Scenario, biasing CV $\psi$

The angle **$\psi$** is a known suboptimal CV for Ala2. It barely distinguishes the two main basins states and is almost orthogonal to the transition state. It is by all means a suboptimal CV.

1.3.1  **Run a single OPES Explore simulation:** 

In analogy with Section 1.2.1, run a standard OPES Explore simulation, this time biasing $\psi$ as a CV, evaluate the FES and visualise the COLVAR trajectory.

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 Are the free energy estimates comparable with the ones biasing $\phi$? What are your impressions? Can one speak of convergence here? </summary>
<br>   
There are very few transitions and the underlying bias is strongly out-of-equilibrium.     
Unsurprisingly, the resulting free energy estimate is very off-target. One cannot extract any meaningful result from this simulation.
</details>  
   
</blockquote>   

One might think that a solution to the problem here is to push the system harder. To test this hypothesis, try increasing the `BARRIER` parameter to 50 kJ/mol or even 100 kJ/mol.

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 Do you see any improvement when increasing `BARRIER`? </summary>
<br>   
One can observe a few more transitions, as one would expected when pushing more, but the bias fluctuations are even more intense. As discussed <a href="#transitions-discussion">before</a>, the free energy estimate is pretty much useless. This is clearly a bad solution that leads again to non converged result.
</details>   
</blockquote>

1.3.2  **Run multiple independent replicas:** 

In analogy with Section 1.2.3, let's try to better assess the convergence quality by running independent simulations with command  
```bash
bash run_multi.sh
```   
To see if longer simulation times allow you to extract a better $\Delta F$ estimate. Try running 5 times longer simulations (50ns) with 
```bash
bash run_multi_50ns.sh
```   

<blockquote> 

[!NOTE]   
<details>
<summary>❓🤔 Can you trust the results and speak of convergence?  </summary>
<br> 
No, the bias is always strongly changing in time making the reweight procedure untrustable, even in the case of longer simulations. The free energy results are not reliable. Averaging over multiple independent simulations confirms the finding from the previous section.
</details>    
</blockquote>


> [!NOTE]
> ❓🤔 Do you have any idea how to improve the quality of the simulations while keeping $\psi$ as the only biased CV?   

---

## Step 2: OneOPES comes to the rescue

OneOPES uses a multi-replica architecture to stimulate an improved sampling. It typically uses 8 replicas (0–7), moving from **convergence-focused** (Replica 0) to **exploration-focused** (Replica 7). Higher exploration-focused replicas have the role to promote transition in phase space by accelerating multiple degrees of freedom, locally with with auxiliary CVs and globally with OPES MultiThermal. The exploratory nature of higher replicas promotes transitions between metastable states that diffuse down the replica ladder through replica exchange. That bias alone may not to be able to build an external quasi-static potential that promotes reversible back and forth transitions. Instead, with the help of higher replicas, replica 0 gains in sampling variety and has a better chance to produce well converged free energy estimates. Note that, for the strategy to be effective, it is paramount that frequent exchanges occur between replicas. The replica exchange algorithm that we use has been implemented in GROMACS by Giovanni Bussi ([MolPhys 2014](http://dx.doi.org/10.1080/00268976.2013.824126)).

### 2.1 A OneOPES Example

Here we propose a OneOPES simulation on Ala2. Replica 0 is analogous to the simulation in Section 1.3 where we put an OPES Explore bias on **$\psi$**.  This is the corresponding plumed.dat script

```plumed
MOLINFO STRUCTURE=input.ala2.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
omega: TORSION ATOMS=@omega-2
ene: ENERGY

d1: DISTANCE ATOMS=2,19
d2: DISTANCE ATOMS=6,17
d3: DISTANCE ATOMS=7,16

OPES_METAD_EXPLORE ...
  LABEL=opes
  ARG=psi
  SIGMA=0.2
  FILE=Kernels.data
  STATE_RFILE=compressed.Kernels
  STATE_WFILE=compressed.Kernels
  PACE=2000
  BARRIER=30
... OPES_METAD_EXPLORE

PRINT STRIDE=500 FILE=COLVAR ARG=opes.bias,phi,psi,ene,omega,d1,d2,d3 FMT=%7.4f

ENDPLUMED
```

There is only one small technical difference betweeen the two scripts: the introduction in OneOPES of an explicit `SIGMA` parameter. `SIGMA` measures the width of the initial basin in CV space that is accessible between PACE depositions. Ideally, `PACE` should be slow enough so that CV diffusion roughly covers all the basins. `SIGMA` is used as a Gaussian width in the OPES probability estimation that leads to the external bias construction process. Too large `SIGMA` values lead to a low resolution external potential that lacks detail greatly hampering its effectiveness. Too small sigmas, can lead to impulsive forces and would take many depositions to fill the phase space. Note that here we are discussing in terms of order of magnitude. Values that are slightly smaller than the initial basing width are typically optimal choices.

While in standard OPES single replica simulations `SIGMA` is estimated by default and it is usual to skip providing it explicitly, in OneOPES we encourage users to estimate it and include it in the script. The estimation can be done with a short unbiased simulation in the initial basin. Even equilibration runs are fine for this purpose. The reason behind this is that its automatic estimation is based on a running average over the dynamics of the CV and this dynamics in replica exchange simulations presents artificially large fluctuations due to coordinate exchanges. The resulting automatic estimate of `SIGMA` would be unnaturally large, lowering the resolution of the OPES external bias and decreasing the overall quality of the simulation.

Replica 1, introduces two additional biases compared to replica 0. It is essential that the differences between replicas are rather limited so that replica exchange acceptances are high.

The first one is a weak OPES Explore bias on an auxiliary CV, in this case the distance between the carbon atoms in the termini of Ala2. This bias has a small `BARRIER` of 3 kJ/mol, comparable to thermal fluctuations, and a `PACE` of 4000, twice slower than the main bias. The second one is an OPES Multithermal bias that lets the system sample temperatures between the thermostat (300 K) and a maximum temperature (here 400 K). `PACE` is set to 500 steps. Instruction `UPDATE_FROM=100` ignores the first 100 ps of simulation for the OPES Multithermal bias construction. This instruction is good practice for large and heterogenous systems such as proteins in membranes, wehere it acts as a mini equilibration before the bias construction and improves its quality. In the Ala2 case, it is not necessary given the tiny size of the system. Here is the corresponding plumed.dat script

```plumed
MOLINFO STRUCTURE=input.ala2.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
omega: TORSION ATOMS=@omega-2
ene: ENERGY

d1: DISTANCE ATOMS=2,19
d2: DISTANCE ATOMS=6,17
d3: DISTANCE ATOMS=7,16

OPES_METAD_EXPLORE ...
  LABEL=opes
  ARG=psi
  SIGMA=0.2
  FILE=Kernels.data
  STATE_RFILE=compressed.Kernels
  STATE_WFILE=compressed.Kernels
  PACE=2000
  BARRIER=30
... OPES_METAD_EXPLORE

OPES_METAD_EXPLORE ...
  LABEL=opes1
  ARG=d1
  SIGMA=0.02
  FILE=Kernels1.data
  STATE_RFILE=compressed1.Kernels
  STATE_WFILE=compressed1.Kernels
  PACE=4000
  BARRIER=3
... OPES_METAD_EXPLORE

ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=400
opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=500 UPDATE_FROM=100

PRINT STRIDE=1000 FILE=COLVAR ARG=opes.bias,phi,psi,ene,omega,d1,d2,d3,opes1.bias,opesX.bias FMT=%7.4f

ENDPLUMED   
```   

In replica 2 we add another weak bias on another distance between different heavy atoms and raise the OPES Multithermal maximum temperature to 500 K. We continue like this until replica 7 that presents in total 3 extra biases on distances and an OPES Multithermal reaching a maximum temperature of 1000 K. This is the corresponding plumed.dat script
   
```plumed
MOLINFO STRUCTURE=input.ala2.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
omega: TORSION ATOMS=@omega-2
ene: ENERGY

d1: DISTANCE ATOMS=2,19
d2: DISTANCE ATOMS=6,17
d3: DISTANCE ATOMS=7,16

OPES_METAD_EXPLORE ...
  LABEL=opes
  ARG=psi
  SIGMA=0.2
  FILE=Kernels.data
  STATE_RFILE=compressed.Kernels
  STATE_WFILE=compressed.Kernels
  PACE=2000
  BARRIER=30
... OPES_METAD_EXPLORE

OPES_METAD_EXPLORE ...
  LABEL=opes1
  ARG=d1
  SIGMA=0.02
  FILE=Kernels1.data
  STATE_RFILE=compressed1.Kernels
  STATE_WFILE=compressed1.Kernels
  PACE=4000
  BARRIER=3
... OPES_METAD_EXPLORE

OPES_METAD_EXPLORE ...
  LABEL=opes2
  ARG=d2
  SIGMA=0.015
  FILE=Kernels2.data
  STATE_RFILE=compressed2.Kernels
  STATE_WFILE=compressed2.Kernels
  PACE=4000
  BARRIER=3
... OPES_METAD_EXPLORE

OPES_METAD_EXPLORE ...
  LABEL=opes3
  ARG=d3
  SIGMA=0.01
  FILE=Kernels3.data
  STATE_RFILE=compressed3.Kernels
  STATE_WFILE=compressed3.Kernels
  PACE=4000
  BARRIER=3
... OPES_METAD_EXPLORE

ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=1000
opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=500 UPDATE_FROM=100

PRINT STRIDE=100 FILE=COLVAR ARG=opes.bias,phi,psi,ene,omega,d1,d2,d3,opesX.bias,opes1.bias,opes2.bias,opes3.bias FMT=%7.4f

ENDPLUMED  
```   

2.1.1  **Run a single OneOPES simulation** 

Enter folder 2.1_OneOPES and run script 

```bash
bash run.sh
```   
to start a OneOPES simulation of Ala2.

The underlying GROMACS command is

```bash
mpirun -n 8 gmx_mpi mdrun -s prd.tpr -plumed plumed.dat -pme cpu -nb cpu -bonded cpu -pin on -nsteps 5000000 -multidir 0 1 2 3 4 5 6 7 -hrex -replex 200 &
```

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 How is the sampling quality and the free energy estimation? </summary>
<br>
There is a large improvement compared to the single replica simulations over $\psi$. The results are in agreement with the reference. As expected an optimal CV is irreplaceable and the best results remain those in single replica over $\phi$.
</details>

</blockquote>

Check the exchange acceptance rate with

```bash
grep -A 23 'Replica exchange statistics' 0/md.log
```

Typically an exchange rate acceptance above 20% is good.

> [!IMPORTANT]
> Pay attention to the key parameter `-replex` that determines every how many steps a replica exchange attempt is performed. We set it to a fast value of 200 steps that corresponds to 1/10 of the OPES Explore main bias `PACE`. We found that this ratio between `PACE` and replex is an optimal choice in OneOPES as it allows a slow construction of the OPES Explore bias in replica 0 while coordinate exchanges occur on a faster schedule. Replica 0 achieves its good phase space sampling in large part because of replica exchange. Its OPES Explore bias over a suboptimal CV can be more of a spectator to the conformational landscape rather than the main driver of transitions. This partially lifts the burden of the main bias in replica 0 to be optimal and prevents the construction of a noisy far out-of-equilibrium OPES Explore potential. The resulting OneOPES external potential is instead built in a more relaxed way that, in turn, allows to perform a good quality reweighting and extract better quality free energies.

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 Try exploring the other extreme of parameter settings, i.e. `--replex 20000` or ten times slower than the main bias `PACE`. What do you observe? </summary>
<br>
The $\Delta F$ estimate significantly worsen. Even if the acceptance rate remains high, the actual number of coordinate exchanges is now too low and insufficient. This way, the suboptimal CVs becomes a limiting factor of the sampling. The free energy estimate is still better than the single replica one over $\psi$, but in OneOPES terms this -replex setting is not an effective choice.
</details>

</blockquote>

Go back to `--replex 200` and try another extreme test. Try disabling altogether the main bias by setting the `PACE` to infinity with
```bash
sed -i "s/PACE=2000/PACE=-1/g" */plumed.dat
```
   
This way, replica 0 becomes unbiased and its sampling relies entirely on transitions from higher replicas. Remember that here we simulate 10 ns, while an unbiased trajectory such as the one in Section 1.1 needed a dynamics of many microseconds to show transitions. 

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 How does OneOPES perform in this case where the main bias is disabled?</summary>
<br>
    
The results are surprisingly good, slightly worse than the standard OneOPES example above, but decent. This gives credit to  the hypothesis that, in this case, sampling is dominated by the higher replicas exploration. <b>Disclaimer:</b> do not read too much into this result. In other more complex examples, this result would not be replicable and the main bias would often play an important role. Ala2 is a toy system with 22 atoms in a vacuum that can be heated up to 1000K without consequences. The same cannot be said for a protein in a membrane, for example.
</details>

</blockquote>

2.1.2  **Run multiple OneOPES simulation** 

In analogy with Section 1.3.2, we provide a better error estimation by running a sequence of 5 independent OneOPES simulations. Enter folder 2.1_OneOPES_multi and run
```bash
bash run_multi.sh
```

<blockquote>

[!NOTE]
<details>
<summary>❓🤔 What are your observations</summary>
<br>
The convergence in the default settings is very robust which is reassuring.
</details>

</blockquote>

### 2.2 Play with OneOPES

The OneOPES strategy is not fixed to an immutable default set of parameters. On the contrary, the user is encouraged to explore and fine-tune their parameter space to find its most effective range for the problem at hand. 

In folder 2.2_Play_with_OneOPES you are given as a starting point the same folder structure as the OneOPES default simulation in Section 2.1 and you are encouraged to change parameters and test the effect of the changes. Given what you learned in the tutorial, use your creativity to change parameters such as `PACE`, `BARRIER`, `replex`, the OPES Multithermal temperature range and evaluate the consequences. Do not forget to always keep an eye to the exchange acceptance rate. 

> [!NOTE]
> ❓🤔 Could you find a set of parameters that further improves the free energy results? 

We can offer some guiding principles on parameter optimisation that come from experience in simulating a variety of systems whose complexity goes beyond Ala2. One common doubt is how aggressive to be in terms `BARRIER`. Our findings is that being overly aggressive does not really pay off. Choosing a `BARRIER` just large enough to trigger back and forth transitions is typically the optimal choice. 

> [!IMPORTANT]
> The `PACE` parameter is another crucial parameter whose optimisation is not trivial and and may seem counterintuitive. Even if it shares a name with the deposition `PACE` in Metadynamics, its optimal usage range is very different. If in Metadynamics the bias potential was built from the ground up, one small Gaussian at a time, it was intuitive that to fill a basin faster, a faster deposition `PACE` interval was needed. On the contrary, in OPES and its variant OPES_EXPLORE, the Gaussians are deposited to estimate the system's probability that in turn is used to build an external potential. A relaxed deposition that lets the system sample well the basins and equilibrate between Gaussian depositions often builds a better quality probability estimate that, in turn, builds a better external potential that brings to convergence faster. `PACE` should be chosen large enough so that this relaxation takes place. As larger systems tend to relax on a slower timescale, they tend to favour slower values of `PACE`. In terms of OneOPES and replica exchange, an optimal `PACE` should be large enough to allow configurations from higher replicas to diffuse to the lower one. When dealing with a new system, it is worthwhile spending time to perform some tests and compare the effect of different `PACE` values. We would recoomend users to start from initial values that have been used in similar systems available on PLUMED NEST and optimise from there. 

Regarding the number of used CVs, typically we would not exceed 3 in the main bias. In larger dimensions, the available phase space to explore and fill with a bias grows to a point that the required simulation time to visit all becomes unreachable. This same reasoning can be applied to other single replica CV-based enhanced sampling methods like Metadynamics and OPES. If one is able to compress the relevant degrees of freedom in an optimised two-dimensional space, this is typically a good solution (see ([JPCL 2025](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.5c02079)) for protein folding). 

Regarding the use of auxiliary CVs, there are a number of solutions. We typically add one extra bias for every higher replica, with a delicate `BARRIER` close to thermal fluctuations and a double `PACE` compared to the one in the main bias. These settings are intended to minimisae changes between neighboring replicas so that replica exchange acceptance is maximised. The temperature range of OPES Multithermal that we would recommend using follows the same principle. Raise the maximum temperature very slowly in the lower replicas and gradually increase the $\Delta T$ between replicas as you reach the higher ones. The key guiding principle is to avoid replica exchange bottlenecks and to always keep an eye on the exchange acceptance rate. When the system is fragile, consider being more delicate or adding walls to prevent it from breaking.

---

## Step 3: More complex applications

In this final part, we present two more complex applications that cannot be run in the time window of a short tutorial. For this reason, we provide their input and output so that you can directly analyse the outcome and in the future you can eventually run the example yourselves.

### 3.1 Trypsin-Benzamidine (Protein-Ligand Binding)

This system is typically used to benchmark protein-ligand binding methods. In an earlier paper ([NatCom 2022](https://www.nature.com/articles/s41467-022-33104-3)), we focused on the role of water, analysing it and encoding it in Machine Learning (ML) CVs such as DeepLDA ([JPCL 2020](https://pubs.acs.org/doi/10.1021/acs.jpclett.0c00535)) and DeepTICA ([PNAS 2021](https://www.pnas.org/doi/10.1073/pnas.2113533118)). In OneOPES, we use a simpler, more intuitive and transferable simulation setup where we do not train any ML CV. We bias as main CVs the funnel projection and radius, as auxiliary CVs we push three different water coordinations and use OPES MultiThermal with a maximum temperature of 370K. In Folder Output you can find the COLVAR files corresponding to 5 independent simulations 250 ns long. In file fes_z_natcom.dat in the parent folder, we include our 1D FES reference over the funnel projection. The reference $\Delta F$ is $26.6 \pm 0.3$ kJ/mol. 

### 3.2 Chignolin (Protein Folding)

The double mutant CLN025 of the Chignolin miniprotein is a well-known 10-residue mini-protein that is widely studied in protein folding method development, especially after the milestone paper by Lindorff-Larsen et al. ([Science 2011](https://www.sciencemag.org/lookup/doi/10.1126/science.1208351)) where for the first time a number of mini-proteins were simulated for 100s of microseconds and their folding landscape was analysed. In that paper, Chignolin was the smallest system. We include file `fes_rmsd_ca_DEShaw.dat` with their reference free energy. 

In OneOPES, we accelerate a rather simple main CV from ([JCP 2018](https://pubs.aip.org/jcp/article/149/19/194113/196500/Folding-a-small-protein-using-harmonic-linear)) that linearly combines 6 contacts via HLDA ([JPCL 2018](http://pubs.acs.org/doi/10.1021/acs.jpclett.8b00733)). As auxiliary CVs we bias the water coordination of the center of the Calphas, the termini distance and the radius of gyration. The masimum temperture reached by OPES Multithermal is 400K. We include an example input folder and as output the COLVAR of 5 independent 400 ns simulations. We observe a slight shift compared to the reference $\Delta F = 3.6 \pm 0.4$ kJ/mol.

To better understand the system and improve the folding OneOPES strategy, in 2025 we developed a set of bioinspired CVs and used them to study again mini-protein folding ([JPCL 2025](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.5c02079)). We accelerate as main CVs a hydrogen-bond focused CV and a sidechain focused CV. As auxiliary CVs we bias the water coordination of 7 heay atoms and we reach a maximum temperature of 420 K in OPES Multithermal. We include the input files in folder `Input2025` and the corresponding COLVAR of 5 independent 1000 ns simulations in `Output2025`. We generated our own reference by running on GROMACS 3 unbiased simulations, 100 us long each. The reference free energy file is in `fes_rmsd_ca_gromacs_3x100us.dat`, the corresponding $\Delta F = 6.9 \pm 0.7$ kJ/mol and the reference trajectory is stored on [Zenodo](https://zenodo.org/records/15583283). Here, the agreement between the OneOPES simulations and the reference is robust.

---

## Troubleshooting

| Issue | Possible Cause & Solution |
| :--- | :--- |
| **No Transitions** | **Possible Cause:** `BARRIER` is too low.<br>**Solution:** Slightly increase the `BARRIER`, for example by about 5 kJ/mol, and run a new simulation. Wait a bit. Remember that larger systems have more degrees of freedom and typically require longer simulation times. Try using a higher `PACE` (slower deposition), especially in large systems. Eventually, consider improving your CV. |
| **The system breaks up or visits high energy irrelevant states** | **Possible Cause:** `BARRIER` is too high<br>**Solution:** Decrease the `BARRIER` until you roughly estimate its minimum value that gives you back and forth transitions. Typically this is the ideal choice for `BARRIER`. <hr> **Possible Cause:** `PACE` is too fast<br>**Solution:** Increase `PACE`. Too fast bias deposition can be troublesome and too aggressive. A slow OPES Explore `PACE` allows a better quality external potential construction typically. <hr> **Possible Cause:** MaxTemp of OPES Multithermal too high<br>**Solution:** Decrease the maximum temperature that OPES Multithermal can reach. Some systems are more fragile than others and one must be delicate. For example biosystems in membranes should not be heated up too much, lower than 350 K, for not risking to break the membrane. |
| **Low Exchange Rates** | **Possible Cause:** Large temperature gaps<br>**Solution:** Use a more gradual temperature ramp across replicas. Be delicate especially on the lower replicas like 1 and 2. <hr> **Possible Cause:** Too aggressive auxiliary CVs bias<br>**Solution:** The same is valid regarding the `BARRIER` parameter on the bias on auxiliary CVs, try to decrease it. |
| **Sparse Low Exchange Rates** | **Possible Cause:** OPES MultiThermal did not generate an effective bias<br>**Solution:** Check the potential energy exploration in the replicas affected by the low exchange rates. If it is not exploring the potential energy range diffusively, but remains stuck around some values for extended portions of time. If that is the case, try changing the OPES Multithermal `PACE` or increase the equilibration period with a longer `UPDATE_FROM`. <hr> **Possible Cause:** Bugs in the code<br>**Solution:** When using PLUMED 2.10 version and GROMACS 2024 there is a bug in how replica exchange is handled that suppresses excahnges and leads to incorrect results. For the moment, we recommend using PLUMED 2.9 and GROMACS 2023 that do not present this bug. |

---

## Current Existing Applications

<div align="center">

<img width="800" height="800" alt="ezgif-49494ac320700ea2" src="https://github.com/obzehn/renders/blob/main/gifs/video_ADRB1.gif?raw=true" />

*Figure 4: OneOPES logo and a sodium ion reaching its binding site in a classA GPCR. Image created and rendered by [Nicola Piasentin](https://github.com/obzehn/renders)*

</div>

*   **Ligand Binding:** [JCTC 2024](https://pubs.acs.org/doi/full/10.1021/acs.jctc.4c01112) and [JPCL 2024](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.4c02352).
*   **Protein Folding:** [JPCL 2025](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.5c02079).
*   **Nucleic acid Dynamics:** [JCIM 2024](https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01166).
*   **Surface Binding:** [Langmuir 2024](https://pubs.acs.org/doi/abs/10.1021/acs.langmuir.4c03301).
*   **GPCR Activation:** [JCTC 2025](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00600), [JPCL 2026](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.5c03834) and [preprint 2026](https://www.biorxiv.org/content/10.64898/2026.06.12.731880).
*   **Ion Channel Modulation:** [bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.02.13.705734).


---

## Future Directions
We hope that this tutorial serves as a **starting point** for a user that is interested in applying OneOPES to the system that they are investigating. OneOPES is a versatile framework that can be combined with a number of CVs, from intuition and physics-based to Machine Learning CVs.

---

## Acknowledgements
V.R. would like to thank Hocine El Khaoudi Enyoury, who tested a preliminary version of the tutorial during his visit to Geneva and provided valuable feedback.
