###
# IMPORTS
###
import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff


###
# PLOTTING FUNCTION FOR SOURCE-SINK MODEL
###
def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )

def plot_FPlines(ax, r, K, d, sigma, p, label, FP_color, N_min, N_max):
    # plot trajectory of FPs
    num = 10000
    motilities = np.linspace(0, 40, num=num)
    N_up_FPs = np.zeros(num)
    N_down_FPs = np.zeros(num)
    for i in range(num):
        nu_Rs = 2 * p * motilities[i]
        nu_Ls = 2 * (1 - p) * motilities[i]
        N_up_FPs[i] = K * (1 - d * (sigma * d + sigma * nu_Rs + nu_Ls) / r / (nu_Ls + sigma * d))
        N_down_FPs[i] = nu_Rs / (nu_Ls + sigma * d) * N_up_FPs[i]
    line = ax.plot(N_up_FPs, N_down_FPs, color=FP_color, label=label)
    add_arrow(line[0], None, 'right', 15, FP_color)
    ax.vlines(K*(1-(d*(2+2*p*(sigma-1)))/(2*r*(1-p))), N_min, N_max, linestyle=':', color=FP_color)
    # add line of r_net = 0
    ax.text(K*(1-(d*(2+2*p*(sigma-1)))/(2*r*(1-p)))+0.01, N_min + 0.02, "$\\overline{r}_{net}(N)=0$", color=FP_color, fontsize=12)
    if 2*d*(1+p*(sigma-1))< r*2*(1-p):
        ax.text(0.03,-0.08,"$\\overline{r}_{net}>0$", color=FP_color, fontsize=12)
    else:
        ax.text(0.03,-0.16, "$\\overline{r}_{net}<0$", color=FP_color, fontsize=12)

###
# MAKE THE PLOT
###
# prepare plot
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4.6), constrained_layout=True)
r = 1
d = [0.25, 0.75]
p = 0.5
sigma = 1
K = 1
FP_colors = ['tab:green', 'tab:orange']
labels = ["source-like FP","sink-like FP"]
N_min = -0.8
N_max = 0.8
# plot nu=infinity axis
if p > 1 / 2:
    ax.plot([(1 - p) / p * N_min, (1 - p) / p * N_max], [N_min, N_max], linestyle=':', color='black')
else:
    ax.plot([N_min, N_max], [p / (1 - p) * N_min, p / (1 - p) * N_max], linestyle=':', color='black')
ax.text(N_max-0.25, N_max-0.3, "$N_{\\downarrow}= N_{\\uparrow}$", fontsize=12)
# plot coordinate axes
ax.vlines(0, N_min, N_max, color='black')
ax.hlines(0, N_min, N_max, color='black')
ax.text(0.05, N_max-0.1, "$N_{\\downarrow}$", fontsize=12)
ax.text(N_max-0.1, -0.1, "$N_{\\uparrow}$", fontsize=12)
# plot FP lines
plot_FPlines(ax,r,K,d[0],sigma,p,labels[0],FP_colors[0], N_min, N_max)
plot_FPlines(ax,r,K,d[1],sigma,p,labels[1],FP_colors[1], N_min, N_max)
# trivial FP
ax.plot(0,0,'ro', label='trivial FP')
# legend
ax.legend()
# set limits and ticks
ax.set_xlim([N_min,N_max])
ax.set_ylim([N_min,N_max])
ax.set_xticks([])
ax.set_yticks([])
plt.show()