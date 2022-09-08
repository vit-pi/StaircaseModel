#####
# IMPORTS
#####

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import OutAnimTot.OutAnimTot_Slow as inpt1
import OutAnimTot.OutAnimTot_Fast as inpt2
#####
# VISUALIZE
#####

# variables
SimLength = inpt1.SimLength   # set this
Interval = int(np.floor(inpt1.StopTime/SimLength))
GenSpaceTot = [inpt1.GenSpaceTot, inpt2.GenSpaceTot]
Times = [inpt1.Times, inpt2.Times]
time_indices = [50, 100, 150, 199]

# prepare plot
fig, axs = plt.subplots(ncols=len(time_indices), nrows=2, sharex=True, sharey=True, figsize=(10,4))
for regime in range(2):
    for i in range(len(time_indices)):
        im = axs[regime, i].imshow(GenSpaceTot[regime][time_indices[i]], vmin=0, vmax=100000, origin='lower', cmap='Greys')
        axs[regime, i].set_xticks([])
        axs[regime, i].set_yticks([])
        maximum = 0
        LD = 0
        for j in range(8):
            if GenSpaceTot[regime][time_indices[i]][j][j] > maximum:
                maximum = GenSpaceTot[regime][time_indices[i]][j][j]
                LD = j
        axs[regime, i].text(0, 6, '$R=$'+str(LD+1), fontsize=12)
        for j in range(7):
            axs[regime, i].hlines(j+1-0.5, j+1-0.5, j+2-0.5, colors='black')
            axs[regime, i].vlines(j+1-0.5, j-0.5, j+1-0.5, colors='black')

# plot labels
axs[0, 0].text(-4, 3, 'low\nmotility\n$(\\nu \\ll \\delta)$', fontsize=12)
axs[1, 0].text(-4, 3, 'high\nmotility\n$(\\nu \\gg \\delta)$', fontsize=12)
for i in range(len(time_indices)):
    axs[0, i].set_title('t = '+str(int(Times[0][time_indices[i]]))+' h')
# finish plot
plt.tight_layout()
plt.show()