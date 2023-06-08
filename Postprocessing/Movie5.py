###
# IMPORTS
###

import NSnapDatabase as sd

###
# PLOT
###
# prepare file
fig_titles = ["Original model", "Stochastic switching", "Density-dependent motility"]
sim_indices = [0, 1, 2]   # type of simulations
version = 1    # version of simulation
max_time = 65
sim_delay = 0.5
fps = 20
staircase = True
wildtype = True

file_paths = ["DeadlyOrig", "DeadlySwitch", "DeadlyDens"]
file_specifiers = ["M10_D3E1", "M20_M1E3_S100_D3E1", "M10_SwM1E3_SwD5Ep4_D3E1"]
zip = "Build/Deadly.zip"
animation_name = "Movie5"

def make_plot(sim_indices):
    simulations = []
    for sim in sim_indices:
        file_path = file_paths[sim]#"Snap_" + file_specifiers[sim] + "_V" + str(version)
        simulations.append(sd.SnapDatabase(file_path, file_specifiers[sim], zip))
    sd.animated_plot(simulations, staircase, wildtype, max_time, sim_delay, ["Slow phenotype", "Fast phenotype"], fig_titles, animation_name)

# make animation
make_plot(sim_indices)