# StaircaseModel
This project provides the code to explore a generalisation of the staircase model for bacterial evolution of antibiotic resistance [1]. This software was created to research the effect of bacterial motility on the evolution of antibiotic resistance [2]. It is instructive to firstly read [2], in order to understand this project. The code consists of two parts: preprocessing (C++) and postprocessing (Python). In the preprocessing, the Next Reaction Method is used to sample the stochastic trajectories of the staircase model. The output can contains various features: snapshots of the statespace at respective times, identified first arrival times corresponding to adaptatation, etc. In postprocessing, this output can be used to create various plots, see [2] for examples. More specific explanation and a demo for each part of the code is given in the respective folders.

# Steps to run the code:
## Preprocessing:
(Tested on Ubuntu 20.04, a supercomputing cluster operated by SLURM with 100 nodes (skylake architecture). Estimated install time: 1 hour. Estimated run time for a single Figure in [2]: 4 days. Estimated run time for a single Movie in [2]: 2 hours.)
1) Prepare a folder od a disk with sufficient free space - usually thousands of output files are created after a submission of a job array into the SLURM.
2) Install GNU C Compiler, for example with a command: module load gcc-9.3.0
3) Compile the files, for example with a command: g++ -std=c++11 -O3 -march=skylake MMainFile.cpp Lattice.cpp Functions.cpp
4) Create a submission bash file for SLURM, job arrays are often needed, for example Bash.sh:

    #!/bin/bash  
    #SBATCH --time=0-18:00:00  
    #SBATCH -n 1  
    #SBATCH --array=0-999  
    FileLocation/CompiledFile.out $SLURM_ARRAY_TASK_ID

5) Put the submission file and compiled file into the prepared folder at FileLocation.
6) Submit the job array to the supercomputer via SLURM, for example with a command: sbatch Bash.sh
7) Once all jobs are finished, zip the folder.

## Postprocessing:
(Tested on Windows 11 Home, Windows 11 Home, Intel(R) Core(TM) i7-8550U CPU, 16 GB RAM, Python and PyCharm installed. Estimated install time: 1 hour. Estimated run time for a single Figure in [2]: 10 minutes. Estimated run time for a single Movie in [2]: 10 minutes.)
1) Place the zip file from preprocessing into the Build folder for postprocessing.
2) Run the postprocessing code for producing the required output from the preprocessing data (such as figure).
3) Your final output is ready.

# References:
[1] R. Hermsen, J. Deris, and T. Hwa, Proceedings of the National Academy of Sciences 109, 10775 (2012), https://www.pnas.org/doi/10.1073/pnas.1117716109.

[2] V. Piskovsky, N. Oliveira, unpublished at the moment

[3] M. A. Gibson and J. Bruck, The journal of physical chemistry A 104, 1876 (2000), https://pubs.acs.org/doi/10.1021/jp993732q.
