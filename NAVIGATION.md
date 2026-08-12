# Tutorial: Mastering Enhanced Sampling with OneOPES

This tutorial was first given at the [Next-Generation Molecular Modeling Summer School 2026](https://github.com/crs4/NGMM2026).

It is a self-contained tutorial that guides the user in performing enhanced sampling simulations, starting from the basics to dealing with suboptimal collective variables and introducing replica exchange.
Users are encouraged to read the main reference paper and to follow the suggested accessory tutorials, from the basics of PLUMED syntax until the OPES tutorial. These tasks are given as inputs in the flow chart below.
A number of OneOPES applications and a tutorial on OneOPES applied to host-guest systems are refereed as outputs in the flow chart.

```mermaid
flowchart TB;
  A[OneOPES paper] ==> B[Instructions];
  C[PLUMED syntax] ==> B;
  D[Statistical errors in MD] ==> B;
  E[Multiple replicas] ==> B;
  F[OPES tutorial] ==> B;
  B ==> G[OneOPES applications];
  B ==> H[OneOPES host-guest tutorial];
  click A "ref1" "The paper where OneOPES was introduced.";
  click B "README.md" "The OneOPES tutorial.";
  click C "ref2" "PLUMED Masterclass 21.1: PLUMED syntax and analysis";
  click D "ref3" "PLUMED Masterclass 21.2: Statistical errors in MD";
  click E "ref4" "PLUMED Masterclass 21.5: Simulations with multiple replicas";
  click F "ref5" "PLUMED Masterclass 22.03: Rethinking Metadynamics: the On-the-fly Probability Enhanced Sampling (OPES) method";
  click G "README.md#Applications" "OneOPES current application list, with their paper reference and their corresponding input files.";
  click H "ref6" "PLUMED Masterclass 24.16: OneOPES Tutorial for Host-Guest Systems with PLUMED";
```

## References
1. Invernizzi, Piaggi, and Parrinello, [Unified Approach to Enhanced Sampling](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034), PRX (2020)
2. Invernizzi and Parrinello, [Exploration vs Convergence Speed in Adaptive-Bias Enhanced Sampling](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00152), JCTC (2022)
3. Rizzi, Aureli, Ansari and Gervasio, [OneOPES, a Combined Enhanced Sampling Method to Rule Them All](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00254), JCTC (2023)
