###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
LD = 2
time_type = 0
min_max = [3e-6, 8e-3]
sim_num = 5   # new simulation number

# Create ploting functions
def res_cost_plot(ax):
    paths = ["AdapRateResCostV" + str(i + 1) + "/" for i in range(sim_num)]
    type = [18 for i in range(sim_num)]
    zip = ["Build/NewBuild_ResCost.zip" for i in range(sim_num)]
    data_set = ad.AdapDatabase(paths, type, zip)
    [min,max,im] = data_set.res_cost_heatmap(LD, time_type, fig, ax, min_max, False, True)
    print(min)
    print(max)


def res_cost_deadly_plot(ax):
    paths = ["AdapRateResCostDeadlyV" + str(i + 1) + "/" for i in range(sim_num)]
    type = [19 for i in range(sim_num)]
    zip = ["Build/NewBuild_ResCost.zip" for i in range(sim_num)]
    data_set = ad.AdapDatabase(paths, type, zip)
    [min,max,im] = data_set.res_cost_heatmap(LD, time_type, fig, ax, min_max, True, False)
    ax.set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
    print(min)
    print(max)

# Make the plot
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 3.5), constrained_layout=True, sharey=True)
res_cost_plot(axs[0])
res_cost_deadly_plot(axs[1])
plt.show()
