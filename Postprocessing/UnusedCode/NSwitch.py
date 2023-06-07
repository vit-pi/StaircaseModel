###
# IMPORTS
###

import matplotlib.pyplot as plt
import NAdapDatabase as ad

###
# PLOT
###
# Choose parameters
LDs = [1, 4, 6]
time_type = 0

# Create data_set access and labels
titles = ["I) high ($L_D/L > \\delta/r$ at t=0)", "II) low ($L_D/L < \\delta/r$ at t=0)"]
folder_num = 5
folder_paths = ["AdapRateSwitchV"+str(i+1)+"/" for i in range(folder_num)]
zip = ["Build/NewBuildV"+str(i+1)+".zip" for i in range(folder_num)]
types = [11 for i in range(folder_num)]
data_set = ad.AdapDatabase(folder_paths, types, zip)
death_params = [1e-1, 3e-1]


# Make the plot
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(4, 6), sharey=True, constrained_layout=True)
plot_list = []
for col in range(2):
    # FIRST ROW
    ax = axs[col]
    # make plot
    for i in range(len(LDs)):
        LD = LDs[i]
        print("Current LD: " + str(LD))
        legend_name = ["$L_D=$"+str(LD)+", simulations", "$L_D=$"+str(LD)+", analytics"]
        if i == 0:
            [mini, maxi, plots] = data_set.switch_adap_plot(death_params[col], LD, time_type, [False, 4, 1e-2, False], [True, True], ax, [2e-5, 2e-2], False, legend_name, True)
            if col == 0:
                plot_list = [plots[0], plots[2], plots[1]]
        else:
            [mini, maxi, plots] = data_set.switch_adap_plot(death_params[col], LD, time_type, [False, 4, 1e-2, False], [False, True], ax, [2e-5, 2e-2], False, legend_name, False)
            if col == 0:
                plot_list.insert(len(plot_list)-1, plots[0])
                plot_list.insert(len(plot_list)-1, plots[1])
    # make x label
    ax.set_ylabel("adaptation $a_2$ (1/h)")
    ax.set_xlabel("switching rate $s$ (1/h)")
    # make title
    ax.set_title(titles[col])
#fig.suptitle("Initial division region:")
# make legend
#for col in range(2):
#    ax = axs[col]
#    box = ax.get_position()
#    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
#labels = ["$L_D=1$, simulations", "$L_D=1$, only", "$L_D=4$, simulations", "$L_D=4$, only", "$L_D=6$, simulations", "$L_D=6$, only", "any $L_D$, only"]
#fig.legend(plot_list, labels, loc="lower center", bbox_to_anchor=(0.5, 0), ncol=4)
#plt.subplots_adjust(bottom=0.25)
# set y label
axs[1].set_xlabel("switching rate $s$ (1/h)")
plt.show()
