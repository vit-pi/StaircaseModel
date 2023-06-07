###
# IMPORTS
###
import pandas as pd
import NSnapDatabase as sd
import matplotlib.pyplot as plt

###
# PLOT
###
# prepare file
fast_slow = False    # choose True for right side of the plot, choose False for left side of the plot
fig_titles = ["Explicit density-dependent motility", "Implicit density-dependent motility"]
sim_indices = [1, 0]   # type of simulations
version = 1    # version of simulation
max_time = 3000
sim_delay = 1
fps = 20
staircase = True
wildtype = True

if fast_slow:
    file_specifiers = ["M3_SwM1E3_SwD5Ep4_D3E1", "M3_SwM1E3_SwD5Ep4_S10_D3E1_SwB1E3"]
    zip = "Build/DensMotSnapFastSlow.zip"
    animation_name = "DensMotSnapFastSlow"
else:
    file_specifiers = ["M1E3_SwM3_SwD5Ep4_D3E1", "M1E3_SwM3_SwD5Ep4_S10_D3E1_SwB1E3"]
    zip = "Build/DensMotSnapSlowFast.zip"
    animation_name = "DensMotSnapSlowFast"

# make plot
simulations = []
for sim in sim_indices:
    file_path = "Snap_" + file_specifiers[sim] + "_V" + str(version)
    simulations.append(sd.SnapDatabase(file_path, file_specifiers[sim], zip))
fig, gen_space_tot = sd.snap_animated_plot(1000, simulations, staircase, wildtype, max_time, sim_delay, ["Low-density phenotype", "High-density phenotype"], fig_titles, animation_name)
df_expl0 = pd.DataFrame(gen_space_tot[0][0])
df_expl0.to_csv("SourceData/Fig8SIf_expl_pop0_fastslow"+str(fast_slow)+".csv", index=False)
df_expl1 = pd.DataFrame(gen_space_tot[0][1])
df_expl1.to_csv("SourceData/Fig8SIf_expl_pop1_fastslow"+str(fast_slow)+".csv", index=False)
df_impl = pd.DataFrame(gen_space_tot[1][0])
df_impl.to_csv("SourceData/Fig8SIf_iml_fastslow"+str(fast_slow)+".csv", index=False)
plt.show()