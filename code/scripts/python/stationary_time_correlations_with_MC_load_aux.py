import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

############ PARAMETERS #############
step = .01
tau_max = 1500
#####################################

fig, ax = plt.subplots(figsize=(10., 6.*1.05))

print("Loading numerical expectations and stds...")
time_corr = np.genfromtxt("data/tc_mc_full", delimiter=',')

# averages = np.load("./data/averages_5040.npy")
# variances = np.load("./data/variances_5040.npy")
# correlators = np.load("./data/correlators_5040.npy")
# pccs = np.load("./data/pccs_5040.npy")

ax.set_xlim([0, tau_max])
ax.set_ylim([-.01, 1.01])
ax.axhline(0, color="black")

ax.set_xlabel(r'$\tau$', fontsize=20)
ax.set_ylabel(r'$\rho_{ij}$ for proteins', fontsize=20)

taus = np.arange(0,step,tau_max)


ax.plot(time_corr[:150000,0], pccs[:, 9:15])
ax.plot(time_corr[:150000,0], pccs[:, 2:9:2])
ax.plot(time_corr[:,0], time_corr[:,1:])


##############################################
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axins = inset_axes(ax, width="50%", height="50%", bbox_to_anchor=(0.35, .25, 1.1, 1.3), bbox_transform=ax.transAxes, loc="lower left")


axins.plot(time_corr[:150000,0], pccs[:, 9:15])
axins.plot(time_corr[:150000,0], pccs[:, 2:9:2])
axins.plot(time_corr[:,0], time_corr[:,1:])

axins.axhline(0, color="black")

x1, x2, y1, y2 = 0, 5, .98, 1.001 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
plt.yticks(visible=True)
plt.xticks(visible=True)

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax, axins, loc1=3, loc2=1, fc="none", ec="black", linewidth=1.5, zorder=500, alpha=.6)
##############################################


plt.tight_layout()

fig.patch.set_facecolor('none')

plt.savefig('/home/oleg/sync/study/ulster/ISRC_presentation_2024/img/time_correlations.png', dpi=300)

plt.show()

fig.clear()
