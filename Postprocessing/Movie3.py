###
# IMPORTS
###

import NSnapDatabase as sd

###
# PLOT
###
# prepare file
fig_titles = ["Low switching", "High switching"]
animation_name = "Movie3"
sim_indices = [0, 1]   # type of simulations
version = 1    # version of simulation
max_time = 2090
sim_delay = 10
fps = 20
staircase = True
wildtype = True

file_specifiers = ["M2_M1E2_S1E3","M2_M1E2_S5",
                   "M1E2", "M1",
                   "M1_C6E1", "M1_C2E1",
                   "M10_C7E1", "M3E1_C3E1",
                   "M2_M1E2_S1E3", "M2_M1E2_S5",
                   "M3_M1E4_S1E3_D3E1", "M3_M1E4_S1E3_D1E1", "M3_M1E4_S5_D3E1", "M3_M1E4_S5_D1E1",
                   "M3_M1_S1E3_D3E1", "M3_M1_S1E3_D1E1", "M3_M1_S5_D3E1", "M3_M1_S5_D1E1",
                   "M1E2_C7E1", "M1_C7E1", "M1E2_C3E1", "M1_C3E1",
                   "M20_D3E1_R1", "M20_D3E1_R2", "M20_D3E1_R3",
                   "M1E2_RC2E1", "M1_RC2E1",
                   "M10_D4E1_R1", "M10_D4E1_R2", "M10_D4E1_R3",
                   "SwM1E2_SwD1Ep3_D1E1", "SwM3_SwD1Ep3_D1E1", "SwM1E2_SwD1Ep4_D1E1", "SwM3_SwD1Ep4_D1E1",
                   "SwM1E2_SwD1Ep3_D3E1", "SwM3_SwD1Ep3_D3E1", "SwM1E2_SwD1Ep4_D3E1", "SwM3_SwD1Ep4_D3E1",
                   "M1E2_M1E4_S1E3_D3E1", "M1E2_M1E4_S1E3_D1E1", "M1E2_M1E4_S5_D3E1", "M1E2_M1E4_S5_D1E1",
                   "SwM3_SwD5Ep4_D1E1", "SwM3_SwD5Ep4_D3E1", "SwM1E2_SwD5Ep4_D1E1", "SwM1E2_SwD5Ep4_D3E1",
                   "M1E2_HGT1E7", "M1E2_HGT1E4", "M1_HGT1E7", "M1_HGT1E4",
                   "M1E2_SD7E1", "M1_SD7E1",]
zip = "Build/SnapshotBuild_SameAdap.zip"  # "Build/SnapshotBuild.zip"

def make_plot(sim_indices):
    simulations = []
    for sim in sim_indices:
        file_path = "Snap_" + file_specifiers[sim] + "_V" + str(version)
        simulations.append(sd.SnapDatabase(file_path, file_specifiers[sim], zip))

    sd.animated_plot(simulations, staircase, wildtype, max_time, sim_delay, ["Slower phenotype", "Faster phenotype"], fig_titles, animation_name)

# make animation
make_plot(sim_indices)