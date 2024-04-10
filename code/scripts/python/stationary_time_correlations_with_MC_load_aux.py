import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

############ PARAMETERS #############
step = .01
tau_max = 1500
#####################################

print("Loading numerical expectations and stds...")
time_corr = np.genfromtxt("data/tc_mc_full", delimiter=',')

averages = np.load("./data/averages_5040.npy")
variances = np.load("./data/variances_5040.npy")
correlators = np.load("./data/correlators_5040.npy")
pccs = np.load("./data/pccs_5040.npy")

plt.xlim([0, tau_max])
plt.ylim([-.01, 1.01])
plt.axhline(0, color="black")
plt.grid()

plt.xlabel(r'$\tau$', fontsize=20)
plt.ylabel('PCC', fontsize=20)

taus = np.arange(0,step,tau_max)
plt.plot(time_corr[:150000,0], pccs[:, 9:15])
plt.plot(time_corr[:150000,0], pccs[:, 2:9:2])
plt.plot(time_corr[:,0], time_corr[:,1:])
plt.show()
