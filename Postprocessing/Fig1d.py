###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
time_type = 1
move = [1e-2, 1]

# Prepare datasets
data_set = ad.AdapDatabase(["AdapRateStandardV1/"], [2], ["Build/NewBuildV1.zip"]) # use same names when zipping preprocessed data

# Plot the theoretical stable wildtype population
fig, axs = plt.subplots(2, 4, figsize=(13,5))
for regime in range(2):
    for R in range(4):
        ax = axs[regime, R]
        data_set.wild_type_plot((move[regime], 1e-1), 1, 0, 0, True, R+1, 0, ax, True, True, True, True)
plt.show()