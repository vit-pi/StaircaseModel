###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
LD = [[1,1],[3,1]]
time_type = [[0,0],[2,0]]
switch_params = [9.5e4, 3.14e4]   # 5 for deadly
move_params = [[1e-4, 3.16], [3.16, 1e-4]]    # low-high and high-low mixed motility
death = 1e-1
run_index = [[0,1],[5,0]]
# Create data_set access
folder_num = 3
folder_paths = ["AdapRateDensSwitchV"+str(i+1)+"/" for i in range(folder_num)]
types = [17 for i in range(folder_num)]
zip = ["Build/NewBuildDensSwitchV"+str(i+1)+".zip" for i in range(folder_num)]
data_set = ad.AdapDatabase(folder_paths, types, zip)

x_labels = [[False,False],[True,True]]
pop_labels = [[True,False],[True,False]]
r_labels = [[False,True],[False,True]]

# Find minimal and maximal adaptation rate
fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(5,4), constrained_layout=True)
for row in range(2):
    for col in range(2):
        col = 1-col     # iterate in columns in reversed order to properly display y-labels
        ax = axs[row, col]
        data_set.wild_type_plot((death, switch_params[row], move_params[col][0], move_params[col][1]), run_index[row][col], 0, 0, True, LD[row][col], time_type[row][col], ax, x_labels[row][col], pop_labels[row][col], r_labels[row][col], False)
plt.show()
