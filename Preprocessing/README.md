# Preprocessing
This code can be used the perform the stochastic simulations of the staircase model by the Next Reaction Method [1]. Note that this is often time-consuming and requires a supercomputer for efficiency. There are three main ingredients for creating a successful main file (starts with M) for a specific set of simulations: Properties structures, Lattice object and Functions to perform simulations.

## Properties (Properties.h)
This contains a C++ structures LatProp (lattice properties) and PopProp (population properties), which carry all important initial data of the staircase model that are required to initialize the Lattice object.

## Lattice (Lattice.h, Lattice.cpp)
This contains the main object Lattice, and two supportive objects: IndexPriorityQueue (data structure required for efficiency of Next Reaction Method) and NewtonRaphson (numerical method for solving the algebraic equations for the initial wild-type profile). The Lattice is initialised by the LatProp and PopProp structures and updated (in the sense of [1]) with the Update method.

## Functions (Functions.h, Functions.cpp)
This contains various functions for running appropriate simulation of the staircase. AdapImage and SnapshotTot are the most important functions. AdapImage saves the snapshots of the lattice at the times of adaptation in separate files and also saves these times in a special file. SnapshotTot saves the snapshots of the lattice at specified regular times.

# Demo
This demo shows how to create Movie 1 in [2]. Figures are created similarly with AdapImage function.

## Instructions:
1) Prepare three folders: Snap_M1E2_V1, Snap_M1_V1, Code. (Can include more versions Snap_M1E2_V2, etc.)
2) Open MSnapshotTot.cpp and check that it is ready to simulate the low motility regime with motility rate PProp[0].Move = 1E-2.
3) Compile with a command: g++ -std=c++11 -O3 -march=skylake MSnapshotTot.cpp Lattice.cpp Functions.cpp.
4) Save the compiled file as LowMot.out in the folder Code.
5) Modify the motility in the MSnapshotTot.cpp file to simulate the high motility regime: PProp[0].Move = 1.
6) Compile with a command: g++ -std=c++11 -O3 -march=skylake MSnapshotTot.cpp Lattice.cpp Functions.cpp.
7) Save the compiled file as HighMot.out in the folder Code.
4) In each folder Snap_M1E2_V1 (resp. Snap_M1_V1) create a file Bash.sh with the following content:

    #!/bin/bash  
    #SBATCH --time=0-18:00:00  
    #SBATCH -n 1  
    Code/LowMot.out (resp. Code/HighMot.out)

6) Navigate to folders Snap_M1E2_V1 (resp. Snap_M1_V1) and submit your job to SLURM with: sbatch Bash.sh.
7) Once both jobs are finished, zip Snap_M1E2_V1, Snap_M1_V1 into a single folder called SnapshotBuild_LowHighMot.zip and move this zipped folder to the folder of Build in Postprocessing.

## Output:
Zipped folder SnapshotBuild_LowHighMot.zip.

## Estimated run time:
5 hours

# References:
[1] M. A. Gibson and J. Bruck, The journal of physical chemistry A 104, 1876 (2000), https://pubs.acs.org/doi/10.1021/jp993732q.

[2] V. Piskovsky, N. Oliveira, unpublished at the moment
