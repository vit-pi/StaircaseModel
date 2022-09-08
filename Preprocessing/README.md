# Preprocessing
This code can be used the perform the stochastic simulations of the staircase model by the Next Reaction Method. Note that this is often time-consuming and requires a supercomputer for efficiency. There are three main ingredients for creating a successful main file (starts with M) for a specific set of simulations: Properties structures, Lattice object and Functions to perform simulations.

## Properties (Properties.h)
This contains a C++ structures LatProp and PopProp, which carry all important initial data of the staircase model that are required to initialize the Lattice object.

## Lattice (Lattice.h, Lattice.cpp)
This contains the main object Lattice, and two supportive objects: IndexPriorityQueue (data structure required for efficiency of Next Reaction Method) and NewtonRaphson (numerical method for solving the algebraic equations for the initial wild-type profile). The Lattice is initialised by the LatProp and PopProp structures and updated 

