###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
LD = 5
time_type = 0
min_max = None
colorbar = True
labels = True

# Create data_set access
data_set = ad.AdapDatabase(["Build/AdapRateChemotaxV1/"], [1], None)
# Make the plot
fig, ax = plt.subplots(figsize=(8, 5))
data_set.chemotax_heatmap(LD, time_type, fig, ax, min_max, colorbar, labels)
plt.show()
