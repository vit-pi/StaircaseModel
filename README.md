# StaircaseModel
This project mathematically explores bacterial evolution of antibiotic resistance [1], with the aim to research the effect of bacterial motility on the evolution of antibiotic resistance [1]. This work builds on the staircase model of bacterial evolution in antibiotic gradients [2]. To understand this project it is instructive to firstly read [1]. The code consists of two parts: preprocessing (C++) and postprocessing (Python). In the preprocessing, the Next Reaction Method is used to sample the stochastic trajectories of the staircase model. The output can contain various features: snapshots of the statespace at respective times, identified first arrival times corresponding to adaptatation, etc. In postprocessing, this output can be used to create various plots, see [1] for examples. More specific explanation and a demo for each part of the code is given in the respective folders. Moreover, a detailed guidance for reproducing all figures in [1] is provided below.

# Steps to run the code:
## Preprocessing:
(Tested on Ubuntu 20.04, a supercomputing cluster operated by SLURM with 100 nodes (skylake architecture). Estimated install time: 1 hour. Estimated run time for a single Figure in [1]: 4 days. Estimated run time for a single Movie in [1]: 2 hours.)
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

# Reproduction Instructions:
The instructions for the choice of parameters to reproduce the Figures of [2] can be found in the Supplementary Text S10 of the paper. Below, you can find a detailed instructions that explain which MMainFile.cpp file is needed in the preprocessing part and which Python file is needed in the postprocessing part to reproduce a particular figure.

| **Figure** | **Preprocessing**        | **Notes for preprocessing**                                             | **Postprocessing** |
|------------|--------------------------|-------------------------------------------------------------------------|--------------------|
| 1C         | MSnapshotTot.cpp         | SimInterval=1; configure with values from Supp. Text S10 (Figure Notes) | Fig1c.py           |
| 1D         | MStandardAdap.cpp        | x                                                                       | Fig1d.py           |
| 2B         | MStandardAdap.cpp        | x                                                                       | Fig2b.py           |
| 2C         | x                        | x                                                                       | Fig2c.py           |
| 2D         | x                        | x                                                                       | Fig2d.py           |
| 3B         | MMotSwitchDeadlyAdap.cpp | x                                                                       | Fig3b.py           |
| 3C         | MMotSwitchDeadlyAdap.cpp | x                                                                       | Fig3c.py           |
| 4B         | MDensSwitchAdap.cpp      | x                                                                       | Fig4b.py           |
| 4C         | MDensSwitchAdap.cpp      | x                                                                       | Fig4c.py           |
| 5B         | x                        | x                                                                       | Fig5b.py           |
| 5C         | x                        | x                                                                       | Fig5c.py           |
| Supp1C     | MStandardAdap.cpp        | x                                                                       | Fig1SIc.py         |
| Supp2      | MResCostAdap.cpp         | Run twice: (1) DeathRate=1E-1, (2) DeathRate=3E-1                       | Fig2SI.py          |
| Supp3      | MMutAdap.cpp             | x                                                                       | Fig3SI.py          |
| Supp4      | MCompBelowStairAdap.cpp  | x                                                                       | Fig4SI.py          |
| Supp5      | MStressDeathAdap.cpp     | x                                                                       | Fig5SI.py          |
| Supp6      | MChemotaxAdap.cpp        | x                                                                       | Fig6SI.py          |
| Supp7A     | MSwitchAdap.cpp          | x                                                                       | Fig7SIa.py         |
| Supp7B     | MMotSwitchAdap.cpp       | x                                                                       | Fig7SIb.py         |
| Supp7C     | x                        | x                                                                       | Fig7SIc.py         |
| Supp8B     | MDensSwitchAdap.cpp      | x                                                                       | Fig4b.py           |
| Supp8C     | MDensSwitchAdap.cpp      | x                                                                       | Fig8SIc.py         |
| Supp8D     | MDensMotSwitchAdap.cpp   | x                                                                       | Fig8SId.py         |
| Supp8E     | MDensMotSwitchAdap.cpp   | x                                                                       | Fig8SIe.py         |
| Supp8F     | MSnapshotTot.cpp         | SimInterval=1; configure with values from Supp. Text S10 (Figure Notes) | Fig8SIf.py         |
| Supp9      | MHGTAdap.cpp             | Run twice: (1) DeathRate=1E-1, (2) DeathRate=3E-1                       | Fig9.py            |
| Movie1     | MSnapshotTot.cpp         | SimInterval=10; configure with values under the Movie in Supp. Text     | Movie1.py          |
| Movie2     | MSnapshotTot.cpp         | SimInterval=10; configure with values under the Movie in Supp. Text     | Movie2.py          |
| Movie3     | MSnapshotTot.cpp         | SimInterval=10; configure with values under the Movie in Supp. Text     | Movie3.py          |
| Movie4     | MSnapshotTot.cpp         | SimInterval=0.1; configure with values under the Movie in Supp. Text    | Movie4.py          |

# References:
[1] V. Piskovsky, N. Oliveira, preprint, https://www.biorxiv.org/content/10.1101/2022.10.21.513270v1

[2] R. Hermsen, J. Deris, and T. Hwa, Proceedings of the National Academy of Sciences 109, 10775 (2012), https://www.pnas.org/doi/10.1073/pnas.1117716109.

[3] M. A. Gibson and J. Bruck, The journal of physical chemistry A 104, 1876 (2000), https://pubs.acs.org/doi/10.1021/jp993732q.
