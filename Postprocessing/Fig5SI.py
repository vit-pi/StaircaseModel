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
LD = 5
time_type = 0
min_max = (1e-5, 1e-2)
colorbar = True
labels = True
sim_num = 5   # new simulation number

# Create ploting functions
def sdeath_plot(ax):
    sdeath_folder_paths = ["AdapRateStressDeathV" + str(i + 1) + "/" for i in range(sim_num)]
    sdeath_type = [10 for i in range(sim_num)]
    zip = ["Build/NewBuildV" + str(i + 1) + ".zip" for i in range(sim_num)]
    data_set_sdeath = ad.AdapDatabase(sdeath_folder_paths, sdeath_type, zip)
    [min, max, adap_rates, stds, im] = data_set_sdeath.stressdeath_heatmap(LD, time_type, fig, ax, min_max, True, labels)
    df_a = pd.DataFrame(adap_rates)
    df_s = pd.DataFrame(stds)
    df_a.to_csv("SourceData/Fig5SI_adap.csv", index=False)
    df_s.to_csv("SourceData/Fig5SI_stds.csv", index=False)

# Make the plot
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 3.5), constrained_layout=True)
sdeath_plot(axs[1])
plt.show()