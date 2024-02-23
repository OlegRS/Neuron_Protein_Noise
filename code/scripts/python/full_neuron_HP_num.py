import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os

# 0) t
# 1) soma__Gene: 1, -0
# 2) soma__mRNA: 1.7959, 1.34011
# 3) soma__Prot: 7004, 1874.81
# 4) d_1__mRNA: 1.29813, 1.13935
# 5) d_1__Prot: 2130.81, 668.535
#   6) s_1_1__Prot: 198.77, 63.7868
#   7) s_1_2__Prot: 198.77, 63.7868
# 8) d_1-1__mRNA: 0.952986, 0.97621
# 9) d_1-1__Prot: 202.58, 146.263
#   10) s_1_1-1__Prot: 18.8974, 14.256
#   11) s_1_1-2__Prot: 188.974, 136.465
# 12) d_1-1__mRNA: 0.952986, 0.97621
# 13) d_1-1__Prot: 892.124, 568.153
#   14) s_1_2-1__Prot: 83.2206, 53.7035
#   15) s_1_2-2__Prot: 83.2206, 53.7035

fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

num_exp = np.genfromtxt("../../bin/exe/expectations", delimiter=',')
num_std = np.genfromtxt("../../bin/exe/variances", delimiter=',')

axs[0].set_ylabel(r'S_1-1', fontsize=20)
axs[0].plot(num_exp[:,0], num_exp[:,6], label="Soma", color='red', alpha=.7)
axs[0].plot(num_exp[:,0], num_exp[:,6] + num_std[:,6], label="Soma", color='red', alpha=.7, linestyle='--')
axs[0].plot(num_exp[:,0], num_exp[:,6] - num_std[:,6], label="Soma", color='red', alpha=.7, linestyle='--')

axs[1].set_ylabel(r'S_1-1_1', fontsize=20)
axs[1].plot(num_exp[:,0], num_exp[:,10], label="Soma", color='red', alpha=.7)
axs[1].plot(num_exp[:,0], num_exp[:,10] + num_std[:,10], label="Soma", color='red', alpha=.7, linestyle='--')
axs[1].plot(num_exp[:,0], num_exp[:,10] - num_std[:,10], label="Soma", color='red', alpha=.7, linestyle='--')

axs[2].set_ylabel(r'S_1-1_2', fontsize=20)
axs[2].plot(num_exp[:,0], num_exp[:,11], label="Soma", color='blue', alpha=.7)
axs[2].plot(num_exp[:,0], num_exp[:,11] + num_std[:,11], label="Soma", color='blue', alpha=.7, linestyle='--')
axs[2].plot(num_exp[:,0], num_exp[:,11] - num_std[:,11], label="Soma", color='blue', alpha=.7, linestyle='--')

axs[3].set_ylabel(r'S_1-2_2', fontsize=20)
axs[3].plot(num_exp[:,0], num_exp[:,15], label="Soma", color='red', alpha=.7)
axs[3].plot(num_exp[:,0], num_exp[:,15] + num_std[:,15], label="Soma", color='red', alpha=.7, linestyle='--')
axs[3].plot(num_exp[:,0], num_exp[:,15] - num_std[:,15], label="Soma", color='red', alpha=.7, linestyle='--')


axs[3].set_xlabel(r'Time in hours', fontsize=18)

for ax in axs:
    # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
    ax.set_xlim([0,200])
    ax.set_ylim(0)
    # ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
