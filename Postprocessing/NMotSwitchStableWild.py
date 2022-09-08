###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
deadly = True
LD = 4
# If deadly region
if deadly:
    switch_params = [1e-3, 5]   # 5 for deadly
    move_params = [[1e-2, 1e-4], [3.16, 1e-4], [3.16, 1]]    # low motility, mixed motility, high motility
    run_index = [[4,0,1],[0,2,1]]
    # Create data_set access
    file_num = 5
    folder_paths = ["AdapRateMotSwitchDeadlyV"+str(i+1)+"/" for i in range(file_num)]
    zip = ["Build/NewBuildV"+str(i+1)+".zip" for i in range(file_num)]
    types = [12 for i in range(file_num)]
# If no deadly region
else:
    switch_params = [1e-3, 1]   # 5 for deadly
    move_params = [[1e-3, 1e-4], [3.16, 1e-4], [3.16, 1]]    # low motility, mixed motility, high motility
    run_index = [[0,5,3],[5,6,3]]
    # Create data_set access
    old_num = 5
    new_num = 5
    old_folder_paths = ["AdapRateMotSwitchV" + str(i + 3) + "/" for i in range(old_num)]
    old_zip = ["Build/BuildV" + str(i + 1) + ".zip" for i in range(old_num)]
    old_types = [0 for i in range(old_num)]
    new_folder_paths = ["AdapRateMotSwitchV" + str(i + 3) + "/" for i in range(new_num)]
    new_zip = ["Build/BuildV" + str(i + 1) + ".zip" for i in range(new_num)]
    new_types = [0 for i in range(new_num)]
    folder_paths = old_folder_paths + new_folder_paths
    types = old_types + new_types
    zip = old_zip + new_zip

x_labels = [[False,False,False],[True,True,True]]
pop_labels = [[True,False,False],[True,False,False]]
r_labels = [[False,False,True],[False,False,True]]

data_set = ad.AdapDatabase(folder_paths, types, zip)

# Find minimal and maximal adaptation rate
fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, sharey=True, figsize=(7,4), constrained_layout=True)
for row in range(2):
    for col in range(3):
        col = 2-col     # iterate in columns in reversed order to properly display y-labels
        ax = axs[row, col]
        data_set.wild_type_plot((move_params[col][0], move_params[col][1], switch_params[row]), run_index[row][col], 0, 1-row, True, LD, 0, ax, x_labels[row][col], pop_labels[row][col], r_labels[row][col], True)
plt.show()
