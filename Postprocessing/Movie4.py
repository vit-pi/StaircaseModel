###
# IMPORTS
###

import matplotlib.pyplot as plt
import NSnapDatabase as sd
import os
import datetime
import numpy as np

###
# PLOT
###
# !!!!! Need to do it in unzipped way, so that no memory error, or Zip just these two!
# !!!!! Also change K to K_S on y_axis
# prepare file
fps = 1
animation_name = "Movie4"
file_specifiers = ["SwM3_SwD1Ep4_D1E1", "SwM3_SwD1Ep4_D3E1"]    # sink and no sink, high motility
versions = [1, 1]    # version of simulation
titles = ["Source-like environment", "Sink-like environment"]
max_time = 25
sim_delay = 0.1
swarm_pop = 0
labels = True

# Create data_set access
zip = "Build/SwarmSinkSensing.zip"
file_paths = ["Snap_" + file_specifiers[i] + "_V" + str(versions[i]) for i in range(2)]
data_set = [sd.SnapDatabase(file_paths[i], file_specifiers[i], zip) for i in range(2)]

# prepare figure properties
wildtype_num = 2
wildtype_width = 5.5
width = wildtype_width*wildtype_num
height = 5
width_ratio = [wildtype_width for i in range(wildtype_num)]
# create figure
fig, axes = plt.subplots(nrows=1, ncols=wildtype_num, figsize=(width, height), gridspec_kw={'width_ratios': width_ratio})
axs = []
for plot in range(2):
    axs.append(axes[plot])
    #ax_twin = axes[plot].twinx()
    axs.append(None)
# prepare folder
date_time = datetime.datetime.now()
folder_specifier = "_"+str(date_time.date().day)+"_"+str(date_time.date().month)+"_"+str(date_time.date().year)+"_"\
                   + str(date_time.time().hour) +"_"+str(date_time.time().minute)
folder_name = "RawMovies\\"+animation_name+folder_specifier
os.mkdir(folder_name)
# draw the figure for each time and take a snapshot
max_times = [min(max_time, data_set[i].stop_time) for i in range(2)]
for plot_index in range(int(max_time/sim_delay)):
    for subplot in range(2):
        sim_index = int(sim_delay/data_set[subplot].sim_interval)*plot_index
        # plot wildtype
        data_set[subplot].wild_type_plot(sim_index, False, False, True, True, False, axs[2*subplot], axs[2*subplot+1], False)
        # change maximum on y_axis to K_S
        K_S = data_set[subplot].p_prop[0].SwarmingDensity
        axs[2*subplot].set_ylim([0, 2*K_S])
        axs[2*subplot].set_yticks([0, K_S/2, K_S, 3*K_S/2, 2*K_S])
        axs[2*subplot].set_yticklabels(["0", "$S/2$", "$S$", "$3S/2$", "$2S$"])
        # read gen_space_tot
        if data_set[subplot].read_gen_space_tot(sim_index) is None:
            [Time, GenSpaceTot] = [np.inf, np.zeros((data_set[subplot].l_prop.PopulationNumber, data_set[subplot].l_prop.GridBound, data_set[subplot].l_prop.GridBound))]
        else:
            [Time, GenSpaceTot] = data_set[subplot].read_gen_space_tot(sim_index)
        # plot title
        if data_set[subplot].l_prop.PopulationNumber == 1 and data_set[subplot].p_prop[0].ConsiderSwarm:
            S = data_set[subplot].S(GenSpaceTot, 0)
            axs[2*subplot].set_title(titles[subplot])#+", S: "+str(S))
        if subplot == 0:
           fig.suptitle('Time'+": {:6.1f}".format(Time)+"h")
    # take a snapshot and clear the canvas
    fig.savefig(folder_name+"\\"+animation_name+"-"+str(plot_index+1)+".png", dpi=300)  # save the figure to file
    for ax in axs:
        if ax is not None:
            for artist in ax.lines + ax.collections + ax.patches:
                artist.remove()
    print("Created frame: "+str(plot_index))

# run ffmpeg in the command prompt
# ffmpeg command such as H:\Můj disk\Summer Research 2021\Staircase model\Simulations\RawMovies\LowMotility_LowHGT_24_3_2022_10_41>ffmpeg -framerate 25 -i LowMotility_LowHGT-%d.png output.mp4
os.system("h:; "
          "cd H:\\Můj disk\\Summer Research 2021\\Staircase model\\Simulations\\"+folder_name+"; "
          "ffmpeg -framerate "+str(fps)+" -i "+animation_name+"-%d.png "+animation_name+".mov")