###
# IMPORTS
###

import NSnapDatabase as sd

###
# PLOT
###
# prepare file
i = 2   # type of simulation
version = 1    # version of simulation
max_time = 10000
sim_delay = 50
fps = 25
staircase = True
wildtype = True

file_specifiers = ["SwM1E2_SwD1Ep3_D1E1", "SwM3_SwD1Ep3_D1E1", "SwM1E2_SwD1Ep4_D1E1", "SwM3_SwD1Ep4_D1E1",
                   "SwM1E2_SwD1Ep3_D3E1", "SwM3_SwD1Ep3_D3E1", "SwM1E2_SwD1Ep4_D3E1", "SwM3_SwD1Ep4_D3E1",
                   "M3_M1E4_S1E3_D3E1", "M3_M1E4_S1E3_D1E1", "M3_M1E4_S5_D3E1", "M3_M1E4_S5_D1E1",
                   "M3_M1_S1E3_D3E1", "M3_M1_S1E3_D1E1", "M3_M1_S5_D3E1", "M3_M1_S5_D1E1",
                   "M1E2_M1E4_S1E3_D3E1", "M1E2_M1E4_S1E3_D1E1", "M1E2_M1E4_S5_D3E1", "M1E2_M1E4_S5_D1E1",
                   "M1E2_C7E1", "M1_C7E1", "M1E2_C3E1", "M1_C3E1",
                   "M1E2", "M1",
                   "SwM3_SwD5Ep4_D1E1", "SwM3_SwD5Ep4_D3E1", "SwM1E2_SwD5Ep4_D1E1", "SwM1E2_SwD5Ep4_D3E1",
                   "M1E2_HGT1E7", "M1E2_HGT1E4", "M1_HGT1E7", "M1_HGT1E4",
                   "M1E2_SD7E1", "M1_SD7E1",]
zip = "Build/SnapshotBuild.zip"

def make_plot(i):
    file_path = "Snap_" + file_specifiers[i] + "_V" + str(version)
    data_set = sd.SnapDatabase(file_path, file_specifiers[i], zip)
    data_set.animated_plot(staircase, wildtype, max_time, sim_delay, fps, ["Slower population", "Faster population"],
                           file_specifiers[i])


# make animation
make_plot(i)