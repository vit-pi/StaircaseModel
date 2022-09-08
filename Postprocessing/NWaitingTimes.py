#####
# IMPORTS
#####

import numpy as np
import matplotlib.pyplot as plt
import NStaircase as sc

#####
# VISUALIZE
#####

# variables
death_rates = [1e-1, 3e-1]
min_nu_exp = -4
max_nu_exp = 1.5
nu_num = 40
LDs = [1, 3, 6]

# prepare arrays
motilities = np.logspace(min_nu_exp, max_nu_exp, num=nu_num)
T1_times = np.zeros((2, 3, nu_num))
D_times = np.zeros((2, 3, nu_num))
maxT1 = 0
maxD = 0
minT1 = 1e6
minD = 1e6
for col in range(2):
    death = death_rates[col]
    for i in range(3):
        LD = LDs[i]
        print('CURRENT LD:'+str(LD))
        for j in range(nu_num):
            # initialize
            nu = motilities[j]
            lat_prop = sc.LatProp()
            pop_prop = [sc.PopProp()]
            pop_prop[0].Move = nu
            lat_prop.LD = LD
            pop_prop[0].DeathRate = death
            stair = sc.Staircase(lat_prop, pop_prop)
            T1 = stair.T1()
            if (not np.isfinite(T1)) or T1 <= 0:
                T1 = 'NaN'
            elif T1>maxT1:
                maxT1 = T1
            elif T1<minT1:
                minT1 = T1
            D = stair.D(1e-2, 1e4)
            if D < 1e-2:
                D = 'NaN'
            elif D>maxD:
                maxD = D
            elif D<minD:
                minD = D
            T1_times[col][i][j] = T1
            D_times[col][i][j] = D

# make plots
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10,3.6), constrained_layout=True, sharex=True, sharey=True)
plot_list = []
for col in range(2):
    ax = axs[col]
    death = death_rates[col]
    # iterate through LD
    for i in range(3):
        LD = LDs[i]
        # plot T1
        mask = np.isfinite(T1_times[col][i])
        pT = ax.plot(motilities[mask], T1_times[col][i][mask])
        plot_list.append(pT[0])
        # plot D
        mask = np.isfinite(D_times[col][i])
        pD = ax.plot(motilities[mask], D_times[col][i][mask], '--', color=pT[0].get_color())
        plot_list.append(pD[0])
    # compute crit_mot and add background
    lat_prop = sc.LatProp()
    pop_prop = [sc.PopProp()]
    lat_prop.LD = 1
    pop_prop[0].DeathRate = death
    stair = sc.Staircase(lat_prop, pop_prop)
    crit_mot = sc.CritMot(10**max_nu_exp, 1e-1, 0, pop_prop, lat_prop)
    if crit_mot is None:
        ax.axvspan(motilities[0]/10, death, facecolor="tab:green", alpha=0.3)
        ax.axvspan(death, motilities[-1]*10, facecolor="tab:orange", alpha=0.3)
        ax.axvline(death, color='black', linestyle="--")
        ax.text(0.4 * death, 1.5 * min(minT1, minD), "$\\delta$", fontsize=12, color='black')
        ax.text(motilities[0]/5, 4 * min(minT1, minD), "Low Motility", fontsize=15, color='green', rotation=90)
        ax.text(2 * motilities[-1], 4 * min(minT1, minD), "High Motility", fontsize=15, color='orange', rotation=90)
    else:
        ax.axvspan(motilities[0]/10, death, facecolor="tab:green", alpha=0.3)
        ax.axvspan(death, crit_mot, facecolor="tab:orange", alpha=0.3)
        ax.axvspan(crit_mot, motilities[-1]*10, facecolor="tab:red", alpha=0.3)
        ax.axvline(death, color='black', linestyle="--")
        ax.text(0.4 * death, 1.5 * min(minT1, minD), "$\\delta$", fontsize=12, color='black')
        ax.axvline(crit_mot, color='black', linestyle="--")
        ax.text(1.2 * crit_mot, 1.5 * min(minT1, minD), "$\\nu_c$", fontsize=12)
        ax.text(motilities[0]/5, 4 * min(minT1, minD), "Low Motility", fontsize=15, color='green', rotation=90)
        ax.text(2 * death, 4 * min(minT1, minD), "High Motility", fontsize=15, color='orange', rotation=90)
        ax.text(2 * motilities[-1], 4 * min(minT1, minD), "Deadly Motility", fontsize=15, color='red', rotation=90)
    # add features
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("motility $\\nu$ (1/h)")
    ax.set_ylim([min(minD, minT1), max(maxT1, maxD)])
    ax.set_xlim([10**min_nu_exp/10, 10**max_nu_exp*10])

# make legend
axs[0].set_ylabel("waiting time (h)")
for col in range(2):
    ax = axs[col]
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
labels = ["$T_1$ time $(L_D=1)$", "$D$ time $(L_D=1)$", "$T_1$ time $(L_D=3)$", "$D$ time $(L_D=3)$", "$T_1$ time $(L_D=6)$", "$D$ time $(L_D=6)$"]
fig.legend(plot_list, labels, loc="lower center", bbox_to_anchor=(0.5, 0), ncol=3)
plt.subplots_adjust(bottom=0.35)
plt.show()