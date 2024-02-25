import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os

# 0) t
# 1) soma__Gene: 1, -0

# 2) soma__mRNA: 1.7959, 1.34011
# 4) d_1__mRNA: 1.29813, 1.13935
# 8) d_1-1__mRNA: 0.952986, 0.97621
# 12) d_1-1__mRNA: 0.952986, 0.97621

# 3) soma__Prot: 7004, 1874.81
# 5) d_1__Prot: 2130.81, 668.535
# 9) d_1-1__Prot: 202.58, 146.263
# 13) d_1-1__Prot: 892.124, 568.153


fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

print("Loading numerical expectations and stds...")
num_exp = np.genfromtxt("data/stationary_expectations", delimiter=',')
num_std = np.genfromtxt("data/stationary_stds", delimiter=',')


################### Gillespie #####################
# 1 - Active genes

# 2 - Soma mRNA
# 4 - d_1 mRNA
# 6 - d_1_1 mRNAs
# 8 - d_1_2 mRNAs

# 3 - Soma proteins
# 5 - d_1 proteins
# 7 - d_1_1 proteins
# 9 - d_1_2 proteins

# 10 - s_1-1 proteins
# 11 - s_1-2 proteins
# 12 - s_1-1_1 proteins
# 13 - s_1-1_2 proteions
# 14 - s_1-2_1 proteins
# 15 - s_1-2_2 proteions

############ PARAMETERS #############
step = .01
n_points = int(7000/step)
# x_lim = n_points*step
sim_run_count_1 = 1 # file_count
sim_run_count_2 = 260 # file_count
n_compartments = 16
#####################################

averages_1 = np.zeros((n_points, n_compartments))
variances_1 = np.zeros((n_points, n_compartments))
for file_id in range(sim_run_count_1):
    print("Loading file " + str(file_id))
    file_name = "../../../data/gillespie/stationary_moments_for_cosyne/SM_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
    averages_1[:, 1:] += data[:, 1:]/sim_run_count_1
    variances_1[:, 1:] += data[:, 1:]**2/sim_run_count_1

averages_2 = np.zeros((n_points, n_compartments))
variances_2 = np.zeros((n_points, n_compartments))
for file_id in range(sim_run_count_2):
    print("Loading file " + str(file_id))
    file_name = "../../../data/gillespie/stationary_moments_for_cosyne/SM_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
    averages_2[:, 1:] += data[:, 1:]/sim_run_count_2
    variances_2[:, 1:] += data[:, 1:]**2/sim_run_count_2
stds_2 = np.sqrt(variances_2 - averages_2**2)

####################### COL_1 ###########################
axs[0,0].axhline(num_exp[1]) # Soma mRNAs
axs[0,0].axhline(num_exp[1] + num_std[1], linestyle='--', color='red')
axs[0,0].axhline(num_exp[1] - num_std[1], linestyle='--', color='red')
axs[0,0].plot(data[:,0], averages_1[:,2], label='Soma', color='red', alpha=.5)

axs[0,0].axhline(num_exp[3]) # d_1 mRNAs
axs[0,0].axhline(num_exp[3] + num_std[3], linestyle='--')
axs[0,0].axhline(num_exp[3] - num_std[3], linestyle='--')
axs[0,0].plot(data[:,0], averages_1[:,4], label='D_1', alpha=.5)

axs[0,0].axhline(num_exp[7]) # d_1_1 mRNAs
axs[0,0].axhline(num_exp[7] + num_std[7], linestyle='--')
axs[0,0].axhline(num_exp[7] - num_std[7], linestyle='--')
axs[0,0].plot(data[:,0], averages_1[:,6], label='D_1-1', alpha=.5)

axs[0,0].axhline(num_exp[11]) # d_1_2 mRNAs
axs[0,0].axhline(num_exp[11] + num_std[11], linestyle='--')
axs[0,0].axhline(num_exp[11] - num_std[11], linestyle='--')
axs[0,0].plot(data[:,0], averages_1[:,8], label='D_1-2', alpha=.5)


axs[1,0].axhline(num_exp[2]) # Soma prot
axs[1,0].axhline(num_exp[2] + num_std[2], linestyle='--', color='red')
axs[1,0].axhline(num_exp[2] - num_std[2], linestyle='--', color='red')
axs[1,0].plot(data[:,0], averages_1[:,3], label='Soma', color='red', alpha=.5)

axs[1,0].axhline(num_exp[4]) # d_1 prot
axs[1,0].axhline(num_exp[4] + num_std[4], linestyle='--')
axs[1,0].axhline(num_exp[4] - num_std[4], linestyle='--')
axs[1,0].plot(data[:,0], averages_1[:,5], label='D_1', alpha=.5)

axs[1,0].axhline(num_exp[8]) # d_1_1 prot
axs[1,0].axhline(num_exp[8] + num_std[8], linestyle='--')
axs[1,0].axhline(num_exp[8] - num_std[8], linestyle='--')
axs[1,0].plot(data[:,0], averages_1[:,7], label='D_1-1', alpha=.5)

axs[1,0].axhline(num_exp[12]) # d_1_2 prot
axs[1,0].axhline(num_exp[12] + num_std[12], linestyle='--')
axs[1,0].axhline(num_exp[12] - num_std[12], linestyle='--')
axs[1,0].plot(data[:,0], averages_1[:,9], label='D_1-2', alpha=.5)


axs[2,0].axhline(num_exp[5]) # s_1_1 prot
axs[2,0].axhline(num_exp[5] + num_std[5], linestyle='--', color='red')
axs[2,0].axhline(num_exp[5] - num_std[5], linestyle='--', color='red')
axs[2,0].plot(data[:,0], averages_1[:,10], label='Soma', color='red', alpha=.5)

axs[2,0].axhline(num_exp[6]) # s_1_2 prot
axs[2,0].axhline(num_exp[6] + num_std[6], linestyle='--', color='red')
axs[2,0].axhline(num_exp[6] - num_std[6], linestyle='--', color='red')
axs[2,0].plot(data[:,0], averages_1[:,11], label='Soma', color='red', alpha=.5)

axs[2,0].axhline(num_exp[9]) # s_1-1_1 prot
axs[2,0].axhline(num_exp[9] + num_std[9], linestyle='--')
axs[2,0].axhline(num_exp[9] - num_std[9], linestyle='--')
axs[2,0].plot(data[:,0], averages_1[:,12], label='D_1', alpha=.5)

axs[2,0].axhline(num_exp[10]) # s_1-1_2 prot
axs[2,0].axhline(num_exp[10] + num_std[10], linestyle='--')
axs[2,0].axhline(num_exp[10] - num_std[10], linestyle='--')
axs[2,0].plot(data[:,0], averages_1[:,13], label='D_1', alpha=.5)

axs[2,0].axhline(num_exp[13]) # s_1-2_1 prot
axs[2,0].axhline(num_exp[13] + num_std[13], linestyle='--')
axs[2,0].axhline(num_exp[13] - num_std[13], linestyle='--')
axs[2,0].plot(data[:,0], averages_1[:,14], label='D_1-1', alpha=.5)

axs[2,0].axhline(num_exp[14]) # s_1-2_2 prot
axs[2,0].axhline(num_exp[14] + num_std[14], linestyle='--')
axs[2,0].axhline(num_exp[14] - num_std[14], linestyle='--')
axs[2,0].plot(data[:,0], averages_1[:,15], label='D_1-2', alpha=.5)



####################### COL_2 ###########################

axs[0,1].axhline(num_exp[1]) # Soma mRNAs
axs[0,1].axhline(num_exp[1] + num_std[1], linestyle='--', color='red')
axs[0,1].axhline(num_exp[1] - num_std[1], linestyle='--', color='red')
axs[0,1].plot(data[:,0], averages_2[:,2], label='Soma', color='red', alpha=.5)

axs[0,1].axhline(num_exp[3]) # d_1 mRNAs
axs[0,1].axhline(num_exp[3] + num_std[3], linestyle='--')
axs[0,1].axhline(num_exp[3] - num_std[3], linestyle='--')
axs[0,1].plot(data[:,0], averages_2[:,4], label='D_1', alpha=.5)

axs[0,1].axhline(num_exp[7]) # d_1_1 mRNAs
axs[0,1].axhline(num_exp[7] + num_std[7], linestyle='--')
axs[0,1].axhline(num_exp[7] - num_std[7], linestyle='--')
axs[0,1].plot(data[:,0], averages_2[:,6], label='D_1-1', alpha=.5)

axs[0,1].axhline(num_exp[11]) # d_1_2 mRNAs
axs[0,1].axhline(num_exp[11] + num_std[11], linestyle='--')
axs[0,1].axhline(num_exp[11] - num_std[11], linestyle='--')
axs[0,1].plot(data[:,0], averages_2[:,8], label='D_1-2', alpha=.5)


axs[1,1].axhline(num_exp[2]) # Soma prot
axs[1,1].axhline(num_exp[2] + num_std[2], linestyle='--', color='red')
axs[1,1].axhline(num_exp[2] - num_std[2], linestyle='--', color='red')
axs[1,1].plot(data[:,0], averages_2[:,3], label='Soma', color='red', alpha=.5)

axs[1,1].axhline(num_exp[4]) # d_1 prot
axs[1,1].axhline(num_exp[4] + num_std[4], linestyle='--')
axs[1,1].axhline(num_exp[4] - num_std[4], linestyle='--')
axs[1,1].plot(data[:,0], averages_2[:,5], label='D_1', alpha=.5)

axs[1,1].axhline(num_exp[8]) # d_1_1 prot
axs[1,1].axhline(num_exp[8] + num_std[8], linestyle='--')
axs[1,1].axhline(num_exp[8] - num_std[8], linestyle='--')
axs[1,1].plot(data[:,0], averages_2[:,7], label='D_1-1', alpha=.5)

axs[1,1].axhline(num_exp[12]) # d_1_2 prot
axs[1,1].axhline(num_exp[12] + num_std[12], linestyle='--')
axs[1,1].axhline(num_exp[12] - num_std[12], linestyle='--')
axs[1,1].plot(data[:,0], averages_2[:,9], label='D_1-2', alpha=.5)


axs[2,1].axhline(num_exp[5]) # s_1_1 prot
axs[2,1].axhline(num_exp[5] + num_std[5], linestyle='--', color='red')
axs[2,1].axhline(num_exp[5] - num_std[5], linestyle='--', color='red')
axs[2,1].plot(data[:,0], averages_2[:,10], label='Soma', color='red', alpha=.5)

axs[2,1].axhline(num_exp[6]) # s_1_2 prot
axs[2,1].axhline(num_exp[6] + num_std[6], linestyle='--', color='red')
axs[2,1].axhline(num_exp[6] - num_std[6], linestyle='--', color='red')
axs[2,1].plot(data[:,0], averages_2[:,11], label='Soma', color='red', alpha=.5)

axs[2,1].axhline(num_exp[9]) # s_1-1_1 prot
axs[2,1].axhline(num_exp[9] + num_std[9], linestyle='--')
axs[2,1].axhline(num_exp[9] - num_std[9], linestyle='--')
axs[2,1].plot(data[:,0], averages_2[:,12], label='D_1', alpha=.5)

axs[2,1].axhline(num_exp[10]) # s_1-1_2 prot
axs[2,1].axhline(num_exp[10] + num_std[10], linestyle='--')
axs[2,1].axhline(num_exp[10] - num_std[10], linestyle='--')
axs[2,1].plot(data[:,0], averages_2[:,13], label='D_1', alpha=.5)

axs[2,1].axhline(num_exp[13]) # s_1-2_1 prot
axs[2,1].axhline(num_exp[13] + num_std[13], linestyle='--')
axs[2,1].axhline(num_exp[13] - num_std[13], linestyle='--')
axs[2,1].plot(data[:,0], averages_2[:,14], label='D_1-1', alpha=.5)

axs[2,1].axhline(num_exp[14]) # s_1-2_2 prot
axs[2,1].axhline(num_exp[14] + num_std[14], linestyle='--')
axs[2,1].axhline(num_exp[14] - num_std[14], linestyle='--')
axs[2,1].plot(data[:,0], averages_2[:,15], label='D_1-2', alpha=.5)


for i in range(2):
    axs[2,i].set_xlabel(r'Time in hours', fontsize=18)
    for ax in axs[:,i]:
        # ax.set_xlim([0,num_exp_1[num_exp_1.shape[0]-1,0]])
        ax.set_xlim([0,3000])
        ax.set_ylim(0)
        # ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
