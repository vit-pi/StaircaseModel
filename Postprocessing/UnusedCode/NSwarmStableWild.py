###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
death = 1e-1
events = [[1,2,3],[4,5,6]]
swarm_move = 1
swarm_density = 6e4
run_index = 0

# If deadly region
# Create data_set access
folder_num = 1
folder_paths = ["AdapRateSwarmDeadlyV"+str(i+1)+"/" for i in range(folder_num)]
zip = ["Build/NewBuildSwarm.zip" for i in range(folder_num)]
types = [14 for i in range(folder_num)]
data_set = ad.AdapDatabase(folder_paths, types, zip)

x_labels = [[False,False,False],[True,True,True]]
pop_labels = [[True,False,False],[True,False,False]]
r_labels = [[False,False,True],[False,False,True]]

# Find minimal and maximal adaptation rate
fig, axs = plt.subplots(ncols=3, nrows=2, sharey=True, figsize=(7,4), constrained_layout=True)
for row in range(2):
    for col in range(3):
        col = 2-col     # iterate in columns in reversed order to properly display y-labels
        ax = axs[row, col]
        data_set.wild_type_plot_swarm((swarm_move,swarm_density,death), run_index, 0, swarm_density, death, True, events[row][col], ax, x_labels[row][col], pop_labels[row][col], r_labels[row][col])
plt.show()
