# Imports
import numpy as np
import matplotlib.pyplot as plt
import copy
import NStaircase as sc

# Parameters
switch_params = [[1e-1, 4e-1, 6.7e-1], [9e-1, 1, 10]]
switch_names = [["$s \leq \\delta$", "$\\delta < s \\lesssim s_c$", "$s = s_c$"], ["$s_c \\lesssim s$", "$s_c < s$", "$s_c \\ll s$"]]
move_num = 3
LD = 1

# Plot heatmap of the deadly region
def PlotDeadlyRegion(switch_param, LD, move_num, ax):
    # prepare x, y axes
    move0_axis = np.zeros((move_num, move_num))
    for m_index in range(move_num):
        move0_axis[m_index] = np.logspace(-3, 3, move_num)
    move1_axis = copy.deepcopy(move0_axis)
    move1_axis = np.transpose(move1_axis)
    # prepare deadly region function
    deadly_indicator = np.zeros((move_num, move_num))
    for i in range(move_num):
        for j in range(move_num):
            move_0 = move0_axis[0][i]
            move_1 = move1_axis[j][0]
            PProp = [sc.PopProp(), sc.PopProp()]
            LProp = sc.LatProp()
            LProp.PopulationNumber = 2
            LProp.LD = LD
            for p_prop in PProp:
                p_prop.DeathRate = 3e-1
                p_prop.SwitchUp = switch_param
                p_prop.SwitchDown = switch_param
            PProp[0].Move = move_0
            PProp[1].Move = move_1
            stair = sc.Staircase(LProp, PProp)
            if stair.N[0] >= 1 or stair.N[8] >= 1:
                deadly_indicator[i][j] = 1
    im = ax.contourf(move0_axis, move1_axis, deadly_indicator, 100, vmin=0, vmax=1, cmap=plt.get_cmap('coolwarm_r'))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(PProp[0].DeathRate, color='black', linestyle="--")
    ax.text(0.3 * PProp[0].DeathRate, 1.4 * 1e-3, "$\\delta$", fontsize=12, color='black')
    ax.axhline(PProp[0].DeathRate, color='black', linestyle="--")
    ax.text(1.4 * 1e-3, 0.3 * PProp[0].DeathRate, "$\\delta$", fontsize=12, color='black')
    return im

# Make the plot
fig, axs = plt.subplots(ncols=3, nrows=2, sharex=True, sharey=True, figsize=(9.8,6))
for row in range(2):
    for col in range(3):
        ax = axs[row, col]
        switch_param = switch_params[row][col]
        ax.set_title(switch_names[row][col])
        im = PlotDeadlyRegion(switch_param, LD, move_num, ax)
# set labels
for col in range(3):
    axs[1, col].set_xlabel("motility $\\nu_1$ (1/h)", fontsize=12)
for row in range(2):
    axs[row, 0].set_ylabel("motility $\\nu_2$ (1/h)", fontsize=12)
# set colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.87, 0.15, 0.015, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax, ticks=[i/10 for i in range(11)])
cbar.set_label("survival probability", fontsize=12, rotation=90)

plt.show()