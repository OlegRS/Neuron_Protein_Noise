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



fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

num_exp = np.genfromtxt("../../bin/exe/expectation_differences_syn_1-1_2", delimiter=',')

axs[0,0].plot(num_exp[:,0], num_exp[:,[1,2,4,8,12]])
axs[1,0].plot(num_exp[:,0], num_exp[:,[10,11]], label=["s_1_1-1", "s_1_1-2"])
axs[2,0].plot(num_exp[:,0], num_exp[:,[6,7,14,15]], label=["s_1-1","s_1-2","s_1_2-1","s_1_2-2"])
axs[3,0].plot(num_exp[:,0], num_exp[:,[3,5,9,13]], label=["soma","d_1","d_1_1","d_1_2"])

for ax in axs[:,0]:
    ax.axvline(1.2e-5*3600, color='black')
    ax.set_xlim([0,2])
    # ax.set_xscale('log')
    # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
    # ax.set_xlim([0,5])
    # ax.set_ylim(0)
    ax.legend(loc=4)

plt.tight_layout()

###################################
num_exp = np.genfromtxt("../../bin/exe/expectation_differences_syn_1_1", delimiter=',')

axs[0,1].plot(num_exp[:,0], num_exp[:,[1,2,4,8,12]])
axs[1,1].plot(num_exp[:,0], num_exp[:,[10,11]], label=["s_1_1-1", "s_1_1-2"])
axs[2,1].plot(num_exp[:,0], num_exp[:,[6,7,14,15]], label=["s_1-1","s_1-2","s_1_2-1","s_1_2-2"])
axs[3,1].plot(num_exp[:,0], num_exp[:,[3,5,9,13]], label=["soma","d_1","d_1_1","d_1_2"])

for ax in axs[:,1]:
    ax.axvline(1.2e-5*3600, color='black')
    ax.set_xlim([0,2])
    # ax.set_xscale('log')
    # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
    # ax.set_xlim([0,5])
    # ax.set_ylim(0)
    ax.legend(loc=4)



################### Gillespie #####################
# 1 - Active genes
# 2 - Soma mRNA
# 3 - Soma proteins
# 4 - d_1 mRNA
# 5 - d_1 proteins
# 6 - d_1_1 mRNAs
# 7 - d_1_1 proteins
# 8 - d_1_2 mRNAs
# 9 - d_1_2 proteins
# 10 - s_1-1 proteins
# 11 - s_1-2 proteins
# 12 - s_1-1_1 proteins
# 13 - s_1-1_2 proteions
# 14 - s_1-2_1 proteins
# 15 - s_1-2_2 proteions

############ PARAMETERS #############
n_points = 9999999
step = .1
x_lim = n_points*step
sim_run_count = 10 # Number of files with Gillespie simulations
mult_sim_run_count = 1 # file_count
n_compartments = 16
#####################################


for file_id in range(mult_sim_run_count):#range(len(os.listdir(main_dir))-1):
    file_name = "../../data/gillespie/heterosynaptic_plasticity_full/HP_syn_12_2_transcription_rate_multiplyer_1" + str(file_id)
    
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[1:n_points+1, 1:]
    averages[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances[:, 1:] += data[:, 1:]**2/mult_sim_run_count







plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

# fig.clear()
