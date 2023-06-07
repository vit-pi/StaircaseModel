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
LD = 4
time_type = 0
min_max = [2e-5, 8e-2] # for standard=False[3e-4, 8e-2]
sim_num = 5   # new simulation number
standard = True
line = False

# Create ploting functions
def hgt_plot(ax):
    paths = ["AdapRateHGTV" + str(i + 1) + "/" for i in range(sim_num)]
    type = [15 for i in range(sim_num)]
    zip = ["Build/NewBuildHGT.zip" for i in range(sim_num)]
    if standard:
        paths = paths+["AdapRateStandardV" + str(i + 1) + "/" for i in range(sim_num)]
        zip = zip+["Build/NewBuildV" + str(i + 1) + ".zip" for i in range(sim_num)]
        type = type+[2 for i in range(sim_num)]
    data_set = ad.AdapDatabase(paths, type, zip)
    [min,max,adap_rates,stds,im] = data_set.hgt_heatmap(LD, time_type, fig, ax, min_max, False, True, line)
    print(min)
    print(max)
    df_a = pd.DataFrame(adap_rates)
    df_s = pd.DataFrame(stds)
    df_a.to_csv("SourceData/Fig9SI_adap_source.csv", index=False)
    df_s.to_csv("SourceData/Fig9SI_stds_source.csv", index=False)



def hgt_deadly_plot(ax):
    paths = ["AdapRateHGTDeadlyV" + str(i + 1) + "/" for i in range(sim_num)]
    type = [16 for i in range(sim_num)]
    zip = ["Build/NewBuildHGT.zip" for i in range(sim_num)]
    if standard:
        paths = paths + ["AdapRateStandardV" + str(i + 1) + "/" for i in range(sim_num)]
        zip = zip + ["Build/NewBuildV" + str(i + 1) + ".zip" for i in range(sim_num)]
        type = type + [2 for i in range(sim_num)]
    data_set = ad.AdapDatabase(paths, type, zip)
    [min,max,adap_rates,stds,im] = data_set.hgt_heatmap(LD, time_type, fig, ax, min_max, True, False, line)
    ax.set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
    print(min)
    print(max)
    df_a = pd.DataFrame(adap_rates)
    df_s = pd.DataFrame(stds)
    df_a.to_csv("SourceData/Fig9SI_adap_sink.csv", index=False)
    df_s.to_csv("SourceData/Fig9SI_stds_sink.csv", index=False)

# Make the plot
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 3.5), constrained_layout=True, sharey=True)
hgt_plot(axs[0])
hgt_deadly_plot(axs[1])
plt.show()