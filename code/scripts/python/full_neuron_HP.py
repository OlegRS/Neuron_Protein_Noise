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

num_exp = np.genfromtxt("../../bin/exe/expectations_s_1_2-2_transcription_rate_multiplyer_1", delimiter=',')
num_std = np.genfromtxt("../../bin/exe/variances_s_1_2-2_transcription_rate_multiplyer_1", delimiter=',')


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
n_points = 1000000
# step = .01
# x_lim = n_points*step
sim_run_count = 10 # Number of files with Gillespie simulations
mult_sim_run_count = 3 # file_count
n_compartments = 16
#####################################


averages = np.zeros((n_points, n_compartments))
variances = np.zeros((n_points, n_compartments))

for file_id in range(mult_sim_run_count):#range(len(os.listdir(main_dir))-1):
    print("Loading file " + str(file_id))
    # file_name = "../../data/gillespie/heterosynaptic_plasticity_full_new/HP_syn_12_2_transcription_rate_multiplyer_1_" + str(file_id+2)
    file_name = "../../data/gillespie/heterosynaptic_plasticity_full_new/HP_syn_12_2_transcription_rate_multiplyer_1_bind_rate_multiplyer_10_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
    averages[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances[:, 1:] += data[:, 1:]**2/mult_sim_run_count



axs[0].plot(data[:,0]-5000+10, averages[:,10], label='Average active genes', color='red', alpha=.5)
axs[1].plot(data[:,0]-5000+10, averages[:,12], label='Average active genes', color='red', alpha=.5)
axs[2].plot(data[:,0]-5000+10, averages[:,13], label='Average active genes', color='red', alpha=.5)
axs[3].plot(data[:,0]-5000+10, averages[:,14], label='Average active genes', color='red', alpha=.5)
axs[3].plot(data[:,0]-5000+10, averages[:,15], label='Average active genes', color='red', alpha=.5)



# axs[0].plot(data[:,0]-5000+10, averages[:,[4,6,8]], color='red', alpha=.7)
# axs[1].plot(data[:,0]-5000+10, averages[:,[3]], color='red', alpha=.7)
# axs[2].plot(data[:,0]-5000+10, averages[:,[5,7,9]], color='red', alpha=.7)
# axs[3].plot(data[:,0]-5000+10, averages[:,[10,11,12,13,14,15]], color='red', alpha=.7)



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


axs[3].set_ylabel(r'S_1-2_1', fontsize=20)
axs[3].plot(num_exp[:,0], num_exp[:,14], label="Soma", color='red', alpha=.7)
axs[3].plot(num_exp[:,0], num_exp[:,14] + num_std[:,14], label="Soma", color='red', alpha=.7, linestyle='--')
axs[3].plot(num_exp[:,0], num_exp[:,14] - num_std[:,14], label="Soma", color='red', alpha=.7, linestyle='--')

axs[3].set_ylabel(r'S_1-2_2', fontsize=20)
axs[3].plot(num_exp[:,0], num_exp[:,15], label="Soma", color='red', alpha=.7)
axs[3].plot(num_exp[:,0], num_exp[:,15] + num_std[:,15], label="Soma", color='red', alpha=.7, linestyle='--')
axs[3].plot(num_exp[:,0], num_exp[:,15] - num_std[:,15], label="Soma", color='red', alpha=.7, linestyle='--')


axs[3].set_xlabel(r'Time in hours', fontsize=18)

for ax in axs:
    # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
    ax.set_xlim([0,500])
    ax.set_ylim(0)
    # ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
