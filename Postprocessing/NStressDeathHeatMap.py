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
data_set = ad.AdapDatabase(["AdapRateStressDeathV1/"], [10], ["Build/NewBuildV1.zip"])
# Make the plot
fig, ax = plt.subplots(figsize=(8, 5))
data_set.stressdeath_heatmap(LD, time_type, fig, ax, min_max, colorbar, labels)
plt.show()