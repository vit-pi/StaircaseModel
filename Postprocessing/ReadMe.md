# Postprocessing
This code can be used to transformed the raw data from preprocessing into meaningful images, see [1] for examples. There are three main objects to be used, each saved within a respective file. These objects are used in rest of the files, which act as the main files for plotting a specific image.
## Staircase (NStaircase.py)
This object can be used for analytical and numerical understanding of the staircase model. It provides an alternative to the time-consuming stochastic simulations obtained in the Preprocessing. This object is vital for the correct functioning of few parts of the AdapDatabase and SnapDatabase objects.
## AdapDatabase (NAdapDatabase.py)
This object can be used for postprocessing of the raw data produced with the AdapImage function in the preprocessing (see Functions.h). When a new type of MMainFile is created in the preprocessing, a new object describing this new experiment must be created. Inspect the end of the current file for examples and follow the pattern. Do not forget to add an instance of this new object into the possible_experiment_types list in the initialization of the AdapDatabase object. The methods of this object are divided into data access and plotting functions. The data access functions can be used to transform the data from the preprocessing into Python objects and to transform the Python objects into meaningful quantities, such as the adaptation rate. The plotting functions are methods which aid the plotting into matplotlib axes.
## SnapDatabase (NAdapDatabase.py)
This object can be used for postprocessing of the raw data produced with the SnapshotTot function in the preprocessing (see Functions.h). The methods of this object are divided into data access and plotting functions. The data access functions can be used to transform the data from the preprocessing into Python objects and to transform these Python objects. The plotting functions are methods which aid the plotting into matplotlib axes.

# References:
