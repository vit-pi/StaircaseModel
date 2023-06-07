###
# IMPORTS
###
import pandas as pd
import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
death = 3e-1
switch_params = [[None, 9.5e4], [1e4, 1e2]]
titles = [["legend", "high $S>K(1-\\delta/r)$"], ["moderate $S \\lesssim K(1-\\delta/r)$", "low $S \\ll K(1-\\delta/r)$"]]
LD = 2
time_type = 0
min_max = None #[1e-7, 3e-2], [4.4543346511333916e-07, 0.005374222986355017] for deadly, [2.779210786895243e-07, 0.008739350751535225] for non-deadly
colorbar = False
labels = False
folder_num = 7

# Create data_set access
folder_paths = ["AdapRateDensSwitchV"+str(i+1)+"/" for i in range(folder_num)]
types = [17 for i in range(folder_num)]
zip = ["Build/NewBuildDensSwitchV"+str(i+1)+".zip" for i in range(folder_num)]
data_set = ad.AdapDatabase(folder_paths, types, zip)

# Find minimal and maximal adaptation rate
fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, figsize=(6.8,5.9))
if min_max is None:
    minimum = 1
    maximum = 1e-4
    for row in range(2):
        for col in range(2):
            ax = axs[row, col]
            switch_param = switch_params[row][col]
            if switch_param is not None:
                [min, max, adap_rates, stds, im] = data_set.dens_switch_heatmap(death, switch_param, LD, time_type, fig, ax, min_max, colorbar, labels)
                if min < minimum:
                    minimum = min
                if max > maximum:
                    maximum = max
                df_a = pd.DataFrame(adap_rates)
                df_s = pd.DataFrame(stds)
                df_a.to_csv("SourceData/Fig4b_adap_row" + str(row) + "_col" + str(col) + ".csv", index=False)
                df_s.to_csv("SourceData/Fig4b_stds_row" + str(row) + "_col" + str(col) + ".csv", index=False)
    min_max = [minimum, maximum]

# Make the plot
for row in range(2):
    for col in range(2):
        ax = axs[row, col]
        switch_param = switch_params[row][col]
        if switch_param is not None:
            ax.set_title(titles[row][col])
            [min, max, adap_rates, stds, im] = data_set.dens_switch_heatmap(death, switch_param, LD, time_type, fig, ax, min_max, colorbar, labels)
            df_a = pd.DataFrame(adap_rates)
            df_s = pd.DataFrame(stds)
            df_a.to_csv("SourceData/Fig4b_adap_row" + str(row) + "_col" + str(col) + ".csv", index=False)
            df_s.to_csv("SourceData/Fig4b_stds_row" + str(row) + "_col" + str(col) + ".csv", index=False)
        else:
            ax.set_xscale("log")
            ax.set_yscale("log")
            delta = 0.1
            ax.axvline(delta, color='black', linestyle="--")
            ax.text(0.3 * delta, 1.4 * 1e-5, "$\\delta$", fontsize=12, color='black')
            ax.axhline(delta, color='black', linestyle="--")
            ax.text(1.4 * 1e-5, 0.3 * delta, "$\\delta$", fontsize=12, color='black')
            ax.text(1e-3, 1e-3, "low\nmotility", horizontalalignment='center')
            ax.text(2, 1e-3, "mixed\nmotility", horizontalalignment='center')
            ax.text(1e-3, 1, "mixed\nmotility", horizontalalignment='center')
            ax.text(2, 1, "high\nmotility", horizontalalignment='center')
# set labels
for col in range(2):
    axs[1, col].set_xlabel("low density motility $\\nu_L$ (1/h)", fontsize=12)
for row in range(2):
    axs[row, 0].set_ylabel("high density motility $\\nu_H$ (1/h)", fontsize=12)
# set colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.87, 0.15, 0.015, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("adaptation rate $a_R$ (1/h)", fontsize=12, rotation=90)
print(min_max)
plt.show()
