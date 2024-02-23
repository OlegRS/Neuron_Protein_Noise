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



fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

num_exp = np.genfromtxt("../../bin/exe/expectation_differences", delimiter=',')

axs[0].plot(num_exp[:,0], num_exp[:,[1,2,4,8,12]])
axs[1].plot(num_exp[:,0], num_exp[:,[10,11]], label=["s_1_1-1", "s_1_1-2"])
axs[2].plot(num_exp[:,0], num_exp[:,[6,7,14,15]], label=["s_1_1","s_1_2","s_1_2-1","s_1_2-2"])
axs[3].plot(num_exp[:,0], num_exp[:,[3,5,9,13]], label=["soma","d_1","d_1-1","d_1-2"])

for ax in axs:
    ax.axvline(1.2e-5*3600, color='black')
    ax.set_xlim([0,50])
    # ax.set_xscale('log')
    # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
    # ax.set_xlim([0,5])
    # ax.set_ylim(0)
    ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
