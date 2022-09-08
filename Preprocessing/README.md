# Preprocessing
This code can be used the perform the stochastic simulations of the staircase model by the Next Reaction Method [1]. Note that this is often time-consuming and requires a supercomputer for efficiency. There are three main ingredients for creating a successful main file (starts with M) for a specific set of simulations: Properties structures, Lattice object and Functions to perform simulations.

## Properties (Properties.h)
This contains a C++ structures LatProp (lattice properties) and PopProp (population properties), which carry all important initial data of the staircase model that are required to initialize the Lattice object.

## Lattice (Lattice.h, Lattice.cpp)
This contains the main object Lattice, and two supportive objects: IndexPriorityQueue (data structure required for efficiency of Next Reaction Method) and NewtonRaphson (numerical method for solving the algebraic equations for the initial wild-type profile). The Lattice is initialised by the LatProp and PopProp structures and updated (in the sense of [1]) with the Update method.

## Functions (Functions.h, Functions.cpp)
This contains various functions for running appropriate simulation of the staircase. AdapImage and SnapshotTot are the most important functions. AdapImage saves the snapshots of the lattice at the times of adaptation in separate files and also saves these times in a special file. SnapshotTot saves the snapshots of the lattice at specified regular times.

# References:
[1] M. A. Gibson and J. Bruck, The journal of physical chemistry A 104, 1876 (2000), https://pubs.acs.org/doi/10.1021/jp993732q.
