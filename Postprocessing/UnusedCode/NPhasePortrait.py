###
# IMPORTS
###
import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff


###
# PLOTTING FUNCTION FOR SOURCE-SINK MODEL
###
def plot_streamlines(ax, r, K, d, sigma, nu, p, FP_color):
    # coordinate arrays
    N_min = min(-2*K*(1-d/r), K*(1-(d*(2+2*p*(sigma-1)))/(2*r*(1-p))), K*(1-(d*(2+2*p*(sigma-1)))/(2*r*(1-p)))*p/(1-p))
    N_max = 2*K*(1-d/r)
    N_source = np.linspace(N_min, N_max, 100)
    N_sink = np.linspace(N_min, N_max, 100)
    # full coorindate arrays
    N_up, N_down = np.meshgrid(N_source, N_sink)
    # plotted vector field
    force_up = r*N_up-r/K*N_up**2-d*N_up-2*nu*p*N_up+2*nu*(1-p)*N_down
    force_down = -sigma*d*N_down+2*nu*p*N_up-2*nu*(1-p)*N_down
    # non-trivial FPs
    nu_R = 2 * p * nu
    nu_L = 2 * (1 - p) * nu
    N_up_FP = K * (1 - d * (sigma * d + sigma * nu_R + nu_L) / r / (nu_L + sigma * d))
    N_down_FP = nu_R / (nu_L + sigma * d) * N_up_FP
    FP_stable = True
    if N_up_FP < 0:
        FP_stable = False
    # plot axes
    ax.vlines(0, N_min, N_max, color='black')
    ax.hlines(0, N_min, N_max, color='black')
    # plot trajectory of FPs
    num = 1000
    motilities = np.linspace(0, 100, num=num)
    N_up_FPs = np.zeros(num)
    N_down_FPs = np.zeros(num)
    for i in range(num):
        nu_Rs = 2 * p * motilities[i]
        nu_Ls = 2 * (1 - p) * motilities[i]
        N_up_FPs[i] = K * (1 - d * (sigma * d + sigma * nu_Rs + nu_Ls) / r / (nu_Ls + sigma * d))
        N_down_FPs[i] = nu_Rs / (nu_Ls + sigma * d) * N_up_FPs[i]
    ax.plot(N_up_FPs, N_down_FPs, linestyle='--', color=FP_color)
    if p>1/2:
        ax.plot([(1-p)/p*N_min, (1-p)/p*N_max], [N_min, N_max], linestyle=':', color='black')
    else:
        ax.plot([N_min, N_max], [p/(1-p)*N_min, p/(1-p)*N_max], linestyle=':', color='black')
    # plot streamplot
    ax.streamplot(N_up, N_down, force_up, force_down, density=(0.5,0.9), color="grey")
    # plot FPs
    if FP_stable:
        ax.plot(0, 0, marker='o', fillstyle='full', markersize=7, color='red')
        ax.plot(N_up_FP, N_down_FP, marker='s', color=FP_color, markersize=7)
    else:
        ax.plot(0, 0, marker='s', markersize=7, color='red')
        ax.plot(N_up_FP, N_down_FP, marker='o', fillstyle='full', color=FP_color, markersize=7)
    # remove ticks and set limits
    ax.set_xlim([N_min, N_max])
    ax.set_ylim([N_min, N_max])
    ax.set_xticks([])
    ax.set_yticks([])


###
# MAKE THE PLOT
###
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(2.5*3,2.5*2), constrained_layout=True)
r = 1
d = [0.25, 0.75]
p = 0.5
sigma = 1
K = 1
nu = [[0, 0.1, 100], [0, 0.1, 5]]
titles = [["$\\nu=0$", "$0<\\nu<\\infty$", "$\\nu \\lesssim \\infty$"],
          ["$\\nu=0$", "$0<\\nu<\\nu_c$", "$\\nu_c<\\nu<\\infty$"]]
FP_colors = ['tab:green', 'tab:orange']
for row in range(2):
    for col in range(3):
        plot_streamlines(axs[row][col],r,K,d[row],sigma,nu[row][col],p,FP_colors[row])
        axs[row][col].set_title(titles[row][col])
plt.show()
