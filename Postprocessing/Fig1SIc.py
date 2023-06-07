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
LDs = [1, 4, 6]
time_type = 1
old_num = 5
new_num = 5
new_new_num = 5

# Create data_set access and labels
titles = ["source-like", "sink-like"]
old_folder_paths = ["AdapRateStandardV"+str(i+1)+"/" for i in range(old_num)]
old_zip = ["Build/BuildV"+str(i+1)+".zip" for i in range(old_num)]
new_folder_paths = ["AdapRateStandardV" + str(i + 1) + "/" for i in range(new_num)]
new_zip = ["Build/NewBuildV" + str(i + 1) + ".zip" for i in range(new_num)]
new_new_folder_paths = ["AdapStandardV" + str(i + 1) + "/" for i in range(new_num)]
new_new_zip = ["Build/NewBuild_StandardV" + str(i + 1) + ".zip" for i in range(new_num)]
folder_paths = old_folder_paths + new_folder_paths
zip = old_zip+new_zip
types = [2 for i in range(old_num+new_num)]
data_set = ad.AdapDatabase(folder_paths, types, zip)
death_params = [1e-1, 3e-1]


# Make the plot
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 5), sharey=True)
plot_list = []
for col in range(2):
    # FIRST ROW
    ax = axs[col]
    # make plot
    for i in range(len(LDs)):
        LD = LDs[i]
        print("Current LD: " + str(LD))
        legend_name = ["$R=$"+str(LD)+", simulations", "$R=$"+str(LD)+", analytics"]
        if i == 0:
            [mini, maxi, adap_rates, stds, plots] = data_set.standard_motility_adap_plot(death_params[col], LD, time_type, [True, 57, 1e-2, True], ax, [2e-5, 2e-2], False, legend_name, True)
            #[mini, maxi, plots] = data_set.standard_motility_adap_plot(death_params[col], LD, time_type, [True, 15, 1e-2, True], ax, [2e-5, 2e-2], False, legend_name, True)
            plot_list = plots
        else:
            [mini, maxi, adap_rates, stds, plots] = data_set.standard_motility_adap_plot(death_params[col], LD, time_type, [True, 57, 1e-2, False], ax, [2e-5, 2e-2], False, legend_name, False)
            #[mini, maxi, plots] = data_set.standard_motility_adap_plot(death_params[col], LD, time_type, [True, 15, 1e-2, False], ax, [2e-5, 2e-2], False, legend_name, False)
            plot_list.insert(len(plot_list)-1, plots[0])
            plot_list.insert(len(plot_list)-1, plots[1])
        df_a = pd.DataFrame(adap_rates)
        df_s = pd.DataFrame(stds)
        df_a.to_csv("SourceData/Fig1SIc_adap_deadly"+str(col)+"_R"+str(LD)+".csv", index=False)
        df_s.to_csv("SourceData/Fig1SIc_stds_deadly"+str(col)+"_R"+str(LD)+".csv", index=False)
    # make x label
    ax.set_xlabel("motility $\\nu$ (1/h)")
    # make title
    ax.set_title(titles[col])
# make legend
for col in range(2):
    ax = axs[col]
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
labels = ["$R=1$, simulations", "$R=1$, analytics", "$R=4$, simulations", "$R=4$, analytics", "$R=6$, simulations", "$R=6$, analytics", "Hermsen, 2012"]
fig.legend(plot_list, labels, loc="lower center", bbox_to_anchor=(0.5, 0), ncol=4)
plt.subplots_adjust(bottom=0.25)
# set y label
axs[0].set_ylabel("adaptation rate $a_R$ (1/h)")
plt.show()