###
# IMPORTS
###

import matplotlib.pyplot as plt
import NStaircase as sc
###
# PLOT
###
# Choose parameters
LD = 2
death_rates = [1e-1, 3e-1]
time_type = 0

# Create data_set access and labels
move_params = [1e-2, 1, 30]
titles = ["$\\nu < \\delta$", "$\\delta < \\nu$", "$\\delta \\ll \\nu$", "$\\nu < \\delta$", "$\\delta < \\nu < \\nu_c$", "$\\nu_c < \\nu$"]

# Prepare plotting data
fig, axs = plt.subplots(ncols=6, nrows=1, constrained_layout=True, figsize=(9, 2))
axs2 = [None for i in range(6)]
# Make plot
for col in range(6):
    ax = axs[col]
    # create staircase object
    move = move_params[col%3]
    l_prop = sc.LatProp()
    p_prop = [sc.PopProp()]
    l_prop.LD = LD
    p_prop[0].Move = move
    p_prop[0].DeathRate = death_rates[col//3]
    stair = sc.Staircase(l_prop, p_prop)
    # make plot
    axs2[col] = ax.twinx()
    ax2 = axs2[col]
    ax2.bar([k + 1 for k in range(stair.LProp.GridBound)], stair.r_net[0], color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    ax.plot([i + 1 for i in range(stair.LProp.GridBound)], stair.N, color='black', label='analytics')
    ax.fill_between([i + 1 for i in range(stair.LProp.GridBound)], stair.N)
    # set ticks
    ax2.set_ylim([-0.31, 0.91])
    ax2.set_yticks([])
    ax.tick_params(axis='y', labelcolor='tab:blue')
    ax.set_xticks([1, LD + 1, stair.LProp.GridBound])
    ax.set_xticklabels(["$1$", "$L_D+1$", "$L$"])
    ax.set_ylim([0, stair.LProp.CarryingCapacity])
    ax.set_xlim([1, stair.LProp.GridBound])
    ax.set_yticks([])
    # set labels
    ax.set_title(titles[col])
    if col == 0:
        K = stair.LProp.CarryingCapacity
        ax.set_ylabel('$N_x$=$\\#$ cells', color='tab:blue')
        ax.set_yticks([0, K / 4, K / 2, 3 * K / 4, K])
        ax.set_yticklabels(["0", "K/4", "K/2", "3K/4", "K"])
    if col == 5:
        ax2.set_ylabel('mutant net birth rate', color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')
        yticks = [-0.3, 0, 0.3, 0.6, 0.9]
        ax2.set_yticks(yticks)
    if col > 0:
        ax.get_shared_x_axes().join(ax, axs[0])
        ax.get_shared_y_axes().join(ax, axs[0])
        ax2.get_shared_y_axes().join(ax2, axs2[0])
# set xlabels
for i in range(6):
    axs[i].set_xlabel("space x")
# set description
plt.show()