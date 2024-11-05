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

# fig, ax = plt.subplots(figsize=(10., 6.*1.05))
# fig, ax = plt.subplots(figsize=(10., 6/1.05))
fig, ax = plt.subplots(figsize=(8/1.2, 6/1.05/1.2))

print("Loading numerical expectations and stds...")
time_corr = np.genfromtxt("data/tc_mc_full", delimiter=',')

averages = np.load("./data/averages_5040.npy")
variances = np.load("./data/variances_5040.npy")
correlators = np.load("./data/correlators_5040.npy")
pccs = np.load("./data/pccs_5040.npy")


ax.set_xlim([0, tau_max])
ax.set_ylim([-.01, 1.01])
ax.axhline(0, color="black")

ax.set_xlabel(r'$\tau\ [$hours$]$', fontsize=18)
# ax.set_ylabel(r'$\rho_{ij}$ for proteins in S_1-2_2', fontsize=18)
ax.set_ylabel(r'Correlations of protein counts', fontsize=17)

taus = np.arange(0,step,tau_max)

# colors = ['#228B22', '#7CFC00', '#4CA64C', '#FFD700', '#FA8072', '#FFBF4C', '#4C4CA6', '#4CFFFF', '#4C4CFF', '#FF4C4C']
colors = ['#228B22', '#7CFC00', '#FFD700', '#FA8072', '#4C4CA6', '#4CFFFF', '#4CA64C', '#FFBF4C', '#4C4CFF', '#FF4C4C']

syn_indices = [7, 8, 10, 11, 13, 14]

for i in range(len(colors)-4):
    ax.plot(time_corr[:,0], time_corr[:,i+1], color=colors[i])
for i in range(6):
    ax.plot(time_corr[:150000,0], pccs[:, 14-i], color=colors[i])
# for i in range(4):
#     ax.plot(time_corr[:150000,0], pccs[:, 2:9:2][:, i], color=colors[9-i])

# eigenvals = [1.7280e+03,4.6166e+01,1.0098e+01,1.2533e+01,1.2407e+01,1.1019e+01,7.9470e+00,7.9083e+00,3.6994e-01,4.2977e-01,4.3560e-03,2.3026e-02,2.3343e-02,1.0535e-01,4.3200e-02,6.1506e-02,8.3333e-02]

eigenvals = [12.4862,10.2006,9.6934,8.2924,9.3113,6.7968,0.004,0.0227,0.0244,0.062,0.1925,0.1112,0.0615,0.0432,0.1667]

for eig in eigenvals:
    plt.axvline(1/eig, color='blue', alpha=.3)

##############################################
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axins = inset_axes(ax, width="50%", height="50%", bbox_to_anchor=(0.35, .25, 1.1, 1.3), bbox_transform=ax.transAxes, loc="lower left")

for i in range(len(colors)-4):
    axins.plot(time_corr[:,0], time_corr[:,i+1], color=colors[i])
for i in range(6):
    axins.plot(time_corr[:150000,0], pccs[:, 14-i], color=colors[i])
# for i in range(4):
#     axins.plot(time_corr[:150000,0], pccs[:, 2:9:2][:, i], color=colors[9-i])


# axins.plot(time_corr[:,0], time_corr[:,1:])

axins.axhline(0, color="black")

x1, x2, y1, y2 = 0, 5, .98, 1.001 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
plt.yticks(visible=True)
plt.xticks(visible=True)

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax, axins, loc1=3, loc2=1, fc="none", ec="black", linewidth=1.5, zorder=500, alpha=.6)
##############################################

for eig in eigenvals:
    plt.axvline(1/eig, color='blue', alpha=.3)

plt.tight_layout()

fig.patch.set_facecolor('none')

# plt.savefig('/home/oleg/sync/study/conferences/ICMNS/poster/img_compressed/time_correlations.png', dpi=300)
plt.savefig('/home/oleg/oleg_windows/olegr/Desktop/study/BBSRC_application/time_correlations.png', dpi=300)

plt.show()

fig.clear()
