# Postprocessing
This code can be used to transformed the raw data from preprocessing into meaningful images, see [1] for examples. There are three main objects to be used, each saved within a respective file. These objects are used in rest of the code files, which act as the main files for plotting a specific image. Notice that the output data from the Preprocessing shall be saved in Build folder in the .zip form, examples of such files are provided.
## Staircase (NStaircase.py)
This object can be used for analytical and numerical understanding of the staircase model. It provides an alternative to the time-consuming stochastic simulations obtained in the Preprocessing. This object is vital for the correct functioning of few parts of the AdapDatabase and SnapDatabase objects. Use the LatProp and PopProp objects to specify the properties of the lattice and respective populations.
## AdapDatabase (NAdapDatabase.py)
This object can be used for postprocessing of the raw data produced with the AdapImage function in the preprocessing (see Functions.h) - primarily used for finding the adaptation times and snapshots at these times. When a new type of MMainFile is created in the preprocessing, a new object describing this new experiment must be created. Inspect the end of the current file for examples and follow the pattern. Do not forget to add an instance of this new object into the possible_experiment_types list in the initialization of the AdapDatabase object. The methods of this object are divided into data access and plotting functions. The data access functions can be used to transform the data from the preprocessing into Python objects and to transform the Python objects into meaningful quantities, such as the adaptation rate. The plotting functions are methods which aid the plotting into matplotlib axes.
## SnapDatabase (NSnapDatabase.py)
This object can be used for postprocessing of the raw data produced with the SnapshotTot function in the preprocessing (see Functions.h) - primarily used to capture snapshots of time-evolution at regular times. The methods of this object are divided into data access and plotting functions. The data access functions can be used to transform the data from the preprocessing into Python objects and to transform these Python objects. The plotting functions are methods which aid the plotting into matplotlib axes.

# Demo
This demo shows how to do Postprocessing of Supplementary Movie 1 in [1]. Preprocessing is covered in the respective folder.

## Instructions:
1) Check that the file SnapshotBuild_LowHighMot.zip, created in Preprocessing, is saved in the folder Build.
2) Create a folder RawMovies.
3) Open the file AnimPlots.py and run the code.
4) Open the folder RawMovies and notice that a folder with the snapshots of the movie were created.
5) Use your favourite software to turn a sequence of movie snapshots into a movie file. (I used ImageJ.)

## Output:
A folder with snapshots of Supplementary Movie 1 in the RawMovies folder. With a use of a third-party software, this can be turned into a movie.

## Estimated run time:
2 minutes

# References:
[1] V. Piskovsky, N. Oliveira, unpublished at the moment
