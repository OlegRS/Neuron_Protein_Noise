import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os

# 1) soma__Gene
# 2) soma__mRNA
# 3) soma__Prot
# 4) d_1__mRNA
# 5) d_1__Prot
# 6) s_1_1__Prot
# 7) s_1_2__Prot
# 8) d_1-1__mRNA
# 9) d_1-1__Prot
# 10) s_1_1-1__Prot
# 11) s_1_1-2__Prot
# 12) d_1-2__mRNA
# 13) d_1-2__Prot
# 14) s_1_2-1__Prot
# 15) s_1_2-2__Prot



fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(10*1.6*1.1, 3.2*1.9*1.8))

num_exp = np.genfromtxt("data/susceptibilities_s_1_1", delimiter=',')

# axs[0,0].plot(num_exp[:,0], num_exp[:,[1,2,4,8,12]])
axs[0,0].set_title('Stationary specific susceptibility to S_1-1 binding rate changes', fontsize=17)
axs[0,0].set_facecolor('lightgray')
axs[0,0].plot(num_exp[:,0], num_exp[:, 6], label= "S_1-1", color='#4CFFFF', linewidth=3, alpha=.7)
axs[1,0].set_facecolor('white')
axs[1,0].plot(num_exp[:,0], num_exp[:,7], label="S_1-2", color='#4C4CA6', linewidth=3, alpha=.7)
axs[2,0].set_facecolor('white')
axs[2,0].plot(num_exp[:,0], num_exp[:,10], label="S_1_1-1", color='#FA8072', linewidth=3, alpha=.7)
axs[2,0].plot(num_exp[:,0], num_exp[:,11], label="S_1_1-2", color='#FFD700', linewidth=3, alpha=.7)
axs[3,0].set_facecolor('white')
axs[3,0].plot(num_exp[:,0], num_exp[:,14], label="S_1_2-1", color='#7CFC00', linewidth=3, alpha=.7)
axs[4,0].set_facecolor('white')
axs[4,0].plot(num_exp[:,0], num_exp[:,15], label="S_1_2-2", color='#228B22', linewidth=3, alpha=.7)
#axs[2,0].plot(num_exp[:,0], num_exp[:,[3,5,9,13]], label=["soma","d_1","d_1_1","d_1_2"])

###################################
num_exp = np.genfromtxt("data/susceptibilities_s_1_2-2", delimiter=',')

# axs[0,1].plot(num_exp[:,0], num_exp[:,[1,2,4,8,12]])
# axs[0,1].plot(num_exp[:,0], num_exp[:,[10,11]], label=["s_1_1-1", "s_1_1-2"])
# axs[1,1].plot(num_exp[:,0], num_exp[:,[6,7,14,15]], label=["s_1-1","s_1-2","s_1_2-1","s_1_2-2"])
# axs[2,1].plot(num_exp[:,0], num_exp[:,[3,5,9,13]], label=["soma","d_1","d_1_1","d_1_2"])
axs[0,1].set_title('Stationary specific susceptibility to S_1_2-2 binding rate changes', fontsize=17)
axs[0,1].plot(num_exp[:,0], num_exp[:, 6], label= "S_1-1", color='#4CFFFF', linewidth=3, alpha=.7)
axs[0,1].set_facecolor('white')
axs[1,1].plot(num_exp[:,0], num_exp[:,7], label="S_1-2", color='#4C4CA6', linewidth=3, alpha=.7)
axs[1,1].set_facecolor('white')
axs[2,1].plot(num_exp[:,0], num_exp[:,10], label="S_1_1-1", color='#FA8072', linewidth=3, alpha=.7)
axs[2,1].plot(num_exp[:,0], num_exp[:,11], label="S_1_1-2", color='#FFD700', linewidth=3, alpha=.7)
axs[2,1].set_facecolor('white')
axs[3,1].plot(num_exp[:,0], num_exp[:,14], label="S_1_2-1", color='#7CFC00', linewidth=3, alpha=.7)
axs[3,1].set_facecolor('white')
axs[4,1].plot(num_exp[:,0], num_exp[:,15], label="S_1_2-2", color='#228B22', linewidth=3, alpha=.7)
axs[4,1].set_facecolor('lightgray')

# axs[0,1].plot(num_exp[:,0], num_exp[:, 6], label= "S_1-1")
# axs[1,1].plot(num_exp[:,0], num_exp[:,7], label="S_1-2")
# axs[2,1].plot(num_exp[:,0], num_exp[:,[10,11]], label=["S_1_1-1", "S_1_1-2"])
# axs[3,1].plot(num_exp[:,0], num_exp[:,14], label="S_1_2-1")
# axs[4,1].plot(num_exp[:,0], num_exp[:,15], label="S_1_2-2")

for i in range(2):
    for ax in axs[:,i]:
        ax.axvline(1.2e-5*3600, color='blue')
        ax.set_xlim([0,1])
        # ax.set_xscale('log')
        # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
        # ax.set_xlim([0,5])
        # ax.set_ylim([-1, 1.7])
        ax.legend(loc='best', fontsize=12)
        ax.axhline(0, color='black')
        ax.tick_params(axis='both', which='major', labelsize=14)
        
        # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for ax in axs[4,:]:
    ax.set_xlabel("Synaptic protein decay rate constant $[$hour$^{-1}]$", fontsize=15)

for ax in axs[:,0]:
    ax.set_ylabel("Susceptibility", fontsize=15)

axs[2,0].set_ylim([-.225, .01])
axs[3,0].set_ylim([-.225, .01])
axs[4,0].set_ylim([-.225, .01])

axs[0,1].set_ylim([-.225, .01])
axs[1,1].set_ylim([-.225, .01])
axs[2,1].set_ylim([-.225, .01])



plt.tight_layout()

fig.patch.set_facecolor('none')
plt.savefig('../../../../../../../conferences/cosyne2024/poster/img/susceptibilities.png', dpi=300)


plt.show()

fig.clear()
