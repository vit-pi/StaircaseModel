#####
# IMPORTS
#####

import numpy as np
import matplotlib.pyplot as plt
import NSnapDatabase as sd

###
# PLOT
###
###
# PLOT
###
# prepare variables
version = 1    # version of simulation
file_specifiers = ["M1E2", "M1"]
time_indices = [250, 500, 750, 995]
# prepare data set
data_set = []
for i in range(2):
    file_path = "Snap_" + file_specifiers[i] + "_V" + str(version)
    zip = "Build/SnapshotBuild_LowHighMot.zip"  # save preprocessing data into this zip folder
    data_set.append(sd.SnapDatabase(file_path, file_specifiers[i], zip))
fig, axs = plt.subplots(ncols=len(time_indices), nrows=2, sharex=True, sharey=True, figsize=(7.5,4))
#make the plot
for i in range(len(time_indices)):
    data_set[0].staircase_plot(time_indices[i], 0, axs[0,i], True)
    data_set[1].staircase_plot(time_indices[i], 0, axs[1, i], True)
# plot labels
axs[0, 0].text(-4, 3, 'low\nmotility\n$(\\nu \\ll \\delta)$', fontsize=12)
axs[1, 0].text(-4, 3, 'high\nmotility\n$(\\nu \\gg \\delta)$', fontsize=12)
# finish plot
plt.tight_layout()
plt.show()