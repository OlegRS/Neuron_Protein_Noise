import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os

fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

num_exp = np.genfromtxt("../../bin/exe/expectations", delimiter=',')
num_std = np.genfromtxt("../../bin/exe/variances", delimiter=',')

axs[0].set_ylabel(r'Active genes', fontsize=20)
axs[0].plot(num_exp[:,0], num_exp[:,1], label="Soma", color='red', alpha=.7)
axs[0].plot(num_exp[:,0], num_exp[:,1] + num_std[:,1], label="Soma", color='red', alpha=.7, linestyle='--')
axs[0].plot(num_exp[:,0], num_exp[:,1] - num_std[:,1], label="Soma", color='red', alpha=.7, linestyle='--')

axs[1].set_ylabel(r'mRNAs', fontsize=20)
axs[1].plot(num_exp[:,0], num_exp[:,2], label="Soma", color='red', alpha=.7)
axs[1].plot(num_exp[:,0], num_exp[:,2] + num_std[:,2], label="Soma", color='red', alpha=.7, linestyle='--')
axs[1].plot(num_exp[:,0], num_exp[:,2] - num_std[:,2], label="Soma", color='red', alpha=.7, linestyle='--')

axs[1].plot(num_exp[:,0], num_exp[:,4], label="Dendrite", color='blue', alpha=.7)
axs[1].plot(num_exp[:,0], num_exp[:,4] + num_std[:,4], label="Dendrite", color='blue', alpha=.7, linestyle='--')
axs[1].plot(num_exp[:,0], num_exp[:,4] - num_std[:,4], label="Dendrite", color='blue', alpha=.7, linestyle='--')

axs[2].set_ylabel(r'Bulk proteins', fontsize=20)
axs[2].plot(num_exp[:,0], num_exp[:,3], label="Soma", color='red', alpha=.7)
axs[2].plot(num_exp[:,0], num_exp[:,3] + num_std[:,3], label="Soma", color='red', alpha=.7, linestyle='--')
axs[2].plot(num_exp[:,0], num_exp[:,3] - num_std[:,3], label="Soma", color='red', alpha=.7, linestyle='--')

axs[2].plot(num_exp[:,0], num_exp[:,5], label="Dendrite", color='blue', alpha=.7)
axs[2].plot(num_exp[:,0], num_exp[:,5] + num_std[:,5], label="Dendrite", color='blue', alpha=.7, linestyle='--')
axs[2].plot(num_exp[:,0], num_exp[:,5] - num_std[:,5], label="Dendrite", color='blue', alpha=.7, linestyle='--')

axs[3].set_ylabel(r'Synaptic proteins', fontsize=18)
axs[3].plot(num_exp[:,0], num_exp[:,6], label="Syn_1-1", color='teal', alpha=.7)
axs[3].plot(num_exp[:,0], num_exp[:,6] + num_std[:,6], label="Syn_1-1", color='teal', alpha=.7, linestyle='--')
axs[3].plot(num_exp[:,0], num_exp[:,6] - num_std[:,6], label="Syn_1-1", color='teal', alpha=.7, linestyle='--')

axs[3].plot(num_exp[:,0], num_exp[:,7], label="Syn_1-2", color='cyan', alpha=.7)
axs[3].plot(num_exp[:,0], num_exp[:,7] + num_std[:,7], label="Syn_1-2", color='cyan', alpha=.7, linestyle='--')
axs[3].plot(num_exp[:,0], num_exp[:,7] - num_std[:,7], label="Syn_1-2", color='cyan', alpha=.7, linestyle='--')

axs[3].set_xlabel(r'Time in hours', fontsize=18)

for ax in axs:
    ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
    ax.set_ylim(0)
    ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
