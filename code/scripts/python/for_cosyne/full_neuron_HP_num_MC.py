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
# 8) d_1-1__mRNA: 0.952986, 0.97621
# 9) d_1-1__Prot: 202.58, 146.263
# 12) d_1-1__mRNA: 0.952986, 0.97621
# 13) d_1-1__Prot: 892.124, 568.153


fig, axs = plt.subplots(nrows=6, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

print("Loading 1st expectations and stds...")
# num_exp_1 = np.genfromtxt("data/expectations_s_1_1_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
# num_std_1 = np.genfromtxt("data/variances_s_1_1_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
# print("Loading 2nd expectations and stds...")
# num_exp_2 = np.genfromtxt("data/expectations_s_1_2-2_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
# num_std_2 = np.genfromtxt("data/variances_s_1_2-2_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
# print("Loading 3rd expectations and stds...")
# num_exp_3 = np.genfromtxt("data/expectations_s_1_2-2_transcription_rate_multiplier_2_bind_r_mult_10", delimiter=',')
# num_std_3 = np.genfromtxt("data/variances_s_1_2-2_transcription_rate_multiplier_2_bind_r_mult_10", delimiter=',')

num_exp_1 = np.genfromtxt("data/2_genes_expectations_s_1_1_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
num_std_1 = np.genfromtxt("data/2_genes_stds_s_1_1_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
print("Loading 2nd expectations and stds...")
num_exp_2 = np.genfromtxt("data/2_genes_expectations_s_1_2-2_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
num_std_2 = np.genfromtxt("data/2_genes_stds_s_1_2-2_transcription_rate_multiplier_1_bind_r_mult_10", delimiter=',')
print("Loading 3rd expectations and stds...")
num_exp_3 = np.genfromtxt("data/2_genes_expectations_s_1_1_transcription_rate_multiplier_2_bind_r_mult_10", delimiter=',')
num_std_3 = np.genfromtxt("data/2_genes_stds_s_1_1_transcription_rate_multiplier_2_bind_r_mult_10", delimiter=',')



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
mult_sim_run_count = 1 # file_count
n_compartments = 16
#####################################

averages_1 = np.zeros((n_points, n_compartments))
variances_1 = np.zeros((n_points, n_compartments))
for file_id in range(mult_sim_run_count):
    print("Loading file " + str(file_id))
    # file_name = "data/HP_syn_1-1_transcription_rate_multiplier_1_bind_rate_multiplier_10_" + str(file_id)
    file_name = "data/HP_syn_1_1_transcription_rate_multiplier_1_bind_rate_multiplier_10_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
    averages_1[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances_1[:, 1:] += data[:, 1:]**2/mult_sim_run_count

averages_2 = np.zeros((n_points, n_compartments))
variances_2 = np.zeros((n_points, n_compartments))
for file_id in range(mult_sim_run_count):
    print("Loading file " + str(file_id))
    # file_name = "data/HP_syn_1-2_2_transcription_rate_multiplier_1_bind_rate_multiplier_10_" + str(file_id)
    file_name = "data/HP_syn_1-2_2_transcription_rate_multiplier_1_bind_rate_multiplier_10_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
    averages_2[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances_2[:, 1:] += data[:, 1:]**2/mult_sim_run_count

averages_3 = np.zeros((n_points, n_compartments))
variances_3 = np.zeros((n_points, n_compartments))
for file_id in range(mult_sim_run_count):
    print("Loading file " + str(file_id))
    # file_name = "data/HP_syn_1-2_2_transcription_rate_multiplier_2_bind_rate_multiplier_10_" + str(file_id)
    # file_name = "data/HP_syn_1-2_2_transcription_rate_multiplier_2_bind_rate_multiplier_10_" + str(file_id)
    file_name = "data/HP_syn_1_1_transcription_rate_multiplier_2_bind_rate_multiplier_10_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
    averages_3[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances_3[:, 1:] += data[:, 1:]**2/mult_sim_run_count


axs[0,0].set_ylabel(r'S_1-1', fontsize=18)
axs[1,0].set_ylabel(r'S_1-2', fontsize=18)
axs[2,0].set_ylabel(r'S_1-1_1', fontsize=18)
axs[3,0].set_ylabel(r'S_1-1_2', fontsize=18)
axs[4,0].set_ylabel(r'S_1-2_1', fontsize=18)
axs[5,0].set_ylabel(r'S_1-2_2', fontsize=18)

######################## 1st column ###########################
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6], label="Expectation", color='red', alpha=.7)
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6] + num_std_1[:,6], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6] - num_std_1[:,6], color='red', alpha=.7, linestyle='--')
axs[0,0].plot(data[:,0]-5000+10, averages_1[:,10], label='Average active genes', color='red', alpha=.5)

axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7], label="Expectation", color='red', alpha=.7)
axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7] + num_std_1[:,7], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7] - num_std_1[:,7], color='red', alpha=.7, linestyle='--')
axs[1,0].plot(data[:,0]-5000+10, averages_1[:,11], label='Average active genes', color='red', alpha=.5)

axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10], label="Expectation", color='red', alpha=.7)
axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10] + num_std_1[:,10], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10] - num_std_1[:,10], color='red', alpha=.7, linestyle='--')
axs[2,0].plot(data[:,0]-5000+10, averages_1[:,12], label='Average active genes', color='red', alpha=.5)

axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11], label="Expectation", color='red', alpha=.7)
axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11] + num_std_1[:,11], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11] - num_std_1[:,11], color='red', alpha=.7, linestyle='--')
axs[3,0].plot(data[:,0]-5000+10, averages_1[:,13], label='Average active genes', color='red', alpha=.5)

axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14], label="Expectation", color='red', alpha=.7)
axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14] + num_std_1[:,14], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14] - num_std_1[:,14], color='red', alpha=.7, linestyle='--')
axs[4,0].plot(data[:,0]-5000+10, averages_1[:,14], label='Average active genes', color='red', alpha=.5)

axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15], label="Expectation", color='red', alpha=.7)
axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15] + num_std_1[:,15], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15] - num_std_1[:,15], color='red', alpha=.7, linestyle='--')
axs[5,0].plot(data[:,0]-5000+10, averages_1[:,15], label='Average active genes', color='red', alpha=.5)

######################## 2nd column ###########################
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6], label="Expectation", color='red', alpha=.7)
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6] + num_std_2[:,6], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6] - num_std_2[:,6], color='red', alpha=.7, linestyle='--')
axs[0,1].plot(data[:,0]-5000+10, averages_2[:,10], label='Average active genes', color='red', alpha=.5)

axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7], label="Expectation", color='red', alpha=.7)
axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7] + num_std_2[:,7], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7] - num_std_2[:,7], color='red', alpha=.7, linestyle='--')
axs[1,1].plot(data[:,0]-5000+10, averages_2[:,11], label='Average active genes', color='red', alpha=.5)

axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10], label="Expectation", color='red', alpha=.7)
axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10] + num_std_2[:,10], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10] - num_std_2[:,10], color='red', alpha=.7, linestyle='--')
axs[2,1].plot(data[:,0]-5000+10, averages_2[:,12], label='Average active genes', color='red', alpha=.5)

axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11], label="Expectation", color='red', alpha=.7)
axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11] + num_std_2[:,11], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11] - num_std_2[:,11], color='red', alpha=.7, linestyle='--')
axs[3,1].plot(data[:,0]-5000+10, averages_2[:,13], label='Average active genes', color='red', alpha=.5)

axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14], label="Expectation", color='red', alpha=.7)
axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14] + num_std_2[:,14], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14] - num_std_2[:,14], color='red', alpha=.7, linestyle='--')
axs[4,1].plot(data[:,0]-5000+10, averages_2[:,14], label='Average active genes', color='red', alpha=.5)

axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15], label="Expectation", color='red', alpha=.7)
axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15] + num_std_2[:,15], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15] - num_std_2[:,15], color='red', alpha=.7, linestyle='--')
axs[5,1].plot(data[:,0]-5000+10, averages_2[:,15], label='Average active genes', color='red', alpha=.5)

######################## 3rd column ###########################
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6], label="Expectation", color='red', alpha=.7)
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6] + num_std_3[:,6], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6] - num_std_3[:,6], color='red', alpha=.7, linestyle='--')
axs[0,2].plot(data[:,0]-5000+10, averages_3[:,10], label='Average active genes', color='red', alpha=.5)

axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7], label="Expectation", color='red', alpha=.7)
axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7] + num_std_3[:,7], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7] - num_std_3[:,7], color='red', alpha=.7, linestyle='--')
axs[1,2].plot(data[:,0]-5000+10, averages_3[:,11], label='Average active genes', color='red', alpha=.5)

axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10], label="Expectation", color='red', alpha=.7)
axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10] + num_std_3[:,10], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10] - num_std_3[:,10], color='red', alpha=.7, linestyle='--')
axs[2,2].plot(data[:,0]-5000+10, averages_3[:,12], label='Average active genes', color='red', alpha=.5)

axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11], label="Expectation", color='red', alpha=.7)
axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11] + num_std_3[:,11], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11] - num_std_3[:,11], color='red', alpha=.7, linestyle='--')
axs[3,2].plot(data[:,0]-5000+10, averages_3[:,13], label='Average active genes', color='red', alpha=.5)

axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14], label="Expectation", color='red', alpha=.7)
axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14] + num_std_3[:,14], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14] - num_std_3[:,14], color='red', alpha=.7, linestyle='--')
axs[4,2].plot(data[:,0]-5000+10, averages_3[:,14], label='Average active genes', color='red', alpha=.5)

axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15], label="Expectation", color='red', alpha=.7)
axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15] + num_std_3[:,15], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15] - num_std_3[:,15], color='red', alpha=.7, linestyle='--')
axs[5,2].plot(data[:,0]-5000+10, averages_3[:,15], label='Average active genes', color='red', alpha=.5)


for i in range(3):
    axs[5,i].set_xlabel(r'Time in hours', fontsize=18)
    for ax in axs[:,i]:
        # ax.set_xlim([0,num_exp_1[num_exp_1.shape[0]-1,0]])
        ax.set_xlim([0,500])
        ax.set_ylim(0)
        # ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
