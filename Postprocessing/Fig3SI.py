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
min_max = [2e-5, 0.2]
sim_num = 5   # new simulation number

paths = ["AdapMutV" + str(i + 1) + "/" for i in range(sim_num)]
type = [20 for i in range(sim_num)]
zip = ["Build/NewBuild_Mut.zip" for i in range(sim_num)]
data_set = ad.AdapDatabase(paths, type, zip)

# Make the plot
fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10, 3.5), constrained_layout=True, sharey=True)
[min,max,adap_rates,stds,im] = data_set.mut_heatmap(0.1, LD, time_type, fig, ax[0], min_max, False, True)
print([min,max])
df_a = pd.DataFrame(adap_rates)
df_s = pd.DataFrame(stds)
df_a.to_csv("SourceData/Fig3SI_adap_source.csv", index=False)
df_s.to_csv("SourceData/Fig3SI_stds_source.csv", index=False)
[min,max,adap_rates,stds,im] = data_set.mut_heatmap(0.3, LD, time_type, fig, ax[1], min_max, True, True)
print([min,max])
ax[1].set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
df_a = pd.DataFrame(adap_rates)
df_s = pd.DataFrame(stds)
df_a.to_csv("SourceData/Fig3SI_adap_sink.csv", index=False)
df_s.to_csv("SourceData/Fig3SI_stds_sink.csv", index=False)
plt.show()