# StaircaseModel
This project provides the code to explore a generalisation of the staircase model for bacterial evolution of antibiotic resistance [1]. This software was created to research the effect of bacterial motility on the evolution of antibiotic resistance [2]. The code consists of two parts: preprocessing (C++) and postprocessing (Python). In the preprocessing, the Next Reaction Method is used to sample the stochastic trajectories of the staircase model. The output can contains various features: snapshots of the statespace at respective times, identified first arrival times corresponding to adaptatation, etc. In postprocessing, this output can be used to create various plots, see [2] for examples.

# Steps to run the code:
1) Preprocessing (works on Linux and requires supercomputing cluster operated by SLURM for efficiency):
  a) Install GNU C Compiler, for example with a command: module load gcc-9.3.0
  b) Compile the files, for example with the code: g++ -std=c++11 -O3 -march=skylake MMainFile.cpp Lattice.cpp Functions.cpp
  b) Create a submission bash file for SLURM. Often job arrays are needed.
2) Postprocessing (works on Windows):

# References:
[1] R. Hermsen, J. Deris, and T. Hwa, Proceedings of the National Academy of Sciences 109, 10775 (2012), https://www.pnas.org/doi/10.1073/pnas.1117716109.

[2] V. Piskovsky, N. Oliveira, unpublished at the moment

[3] M. A. Gibson and J. Bruck, The journal of physical chemistry A 104, 1876 (2000), https://pubs.acs.org/doi/10.1021/jp993732q.
