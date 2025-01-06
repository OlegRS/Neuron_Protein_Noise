import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os

main_dir = "../../data/gillespie/heterosynaptic_plasticity"

file_count = len(os.listdir(main_dir))


############ PARAMETERS #############
n_points = 99999
step = .1
x_lim = n_points*step
sim_run_count = 10 # Number of files with Gillespie simulations
mult_sim_run_count = 1# file_count
dim = 8
#####################################

fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

num_exp = np.genfromtxt("../../data/gillespie/heterosynaptic_plasticity/expectations", delimiter=',')[:n_points, :]

# 0: t

# 0: time
# 1: Soma_AG
# 2: Soma_mRNA
# 3: Soma_Prot
# 4: d_1_mRNA
# 5: d_1_Prot
# 6: s_1_1_Prot
# 7: s_1_1_Prot

averages = np.zeros((n_points, dim))
variances = np.zeros((n_points, dim))

for file_id in range(mult_sim_run_count):#range(len(os.listdir(main_dir))-1):
    file_name = "../../data/gillespie/heterosynaptic_plasticity/HM_" + str(file_id)
    print(file_name)
    data = np.vstack((np.genfromtxt(file_name, delimiter=',')[1:50000, 1:], np.genfromtxt(file_name, delimiter=',')[50001:n_points+2, 1:]))
    
    print("data: \n", data[50000])
    # for i in range(data.shape[0]):
    #     if data[i,0] == "t" :
    #         print("i = ", i)

    averages[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances[:, 1:] += data[:, 1:]**2/mult_sim_run_count
    
#   axs[2].plot(data[:,0], data[:,15], label="Spine_1-2_2", alpha=.5, color='pink')

# for i in range(1,6):
#     for file_id in range(len(os.listdir("../../data/gillespie/single_Y_fork_aux" + str(i)))):
#         file_name = "../../data/gillespie/single_Y_fork_aux" + str(i) + "/single_Y_fork" + str(file_id)
#         print(file_name)
#         data = np.genfromtxt(file_name, delimiter=',')[1:n_points+1, 1:]
#         averages[:, 1:] += data[:, 1:]/mult_sim_run_count
#         variances[:, 1:] += data[:, 1:]**2/mult_sim_run_count


averages[:, 0] = data[:, 0]
variances[:, 0] = data[:, 0]
variances[:, 1:] = np.sqrt(variances[:, 1:] - averages[:, 1:]**2)


# axs[0].plot(num[:,0], num[:,1], color='red', linewidth=2, zorder=1e4)
# axs[0].plot(num[:,0], num[:,1] + num[:,2], linestyle='--', color='red', linewidth=2)
# axs[0].plot(num[:,0], num[:,1] - num[:,2], linestyle='--', color='red', linewidth=2)


axs[0].set_ylabel(r'Active genes', fontsize=20)


axs[0].plot(averages[:,0], averages[:,1], label='Average active genes', color='red', alpha=.5)
axs[0].plot(averages[:,0], averages[:,1] + variances[:,1], linestyle='--', color='red', linewidth=2)
axs[0].plot(averages[:,0], averages[:,1] - variances[:,1], linestyle='--', color='red', linewidth=2)


axs[1].plot(averages[:,0], averages[:,2], label='Soma', color='red', alpha=.5)
axs[1].plot(averages[:,0], averages[:,2] + variances[:,2], linestyle='--', color='red', linewidth=2)
axs[1].plot(averages[:,0], averages[:,2] - variances[:,2], linestyle='--', color='red', linewidth=2)
axs[1].plot(averages[:,0], averages[:,4], label="Dend_1", color='blue', alpha=.5)
axs[1].plot(averages[:,0], averages[:,4] + variances[:,4], linestyle='--', color='blue', linewidth=2)
axs[1].plot(averages[:,0], averages[:,4] - variances[:,4], linestyle='--', color='blue', linewidth=2)
axs[1].plot(averages[:,0], averages[:,6], label="Dend_1-1", color='orange', alpha=.5)
axs[1].plot(averages[:,0], averages[:,6] + variances[:,6], linestyle='--', color='orange', linewidth=2)
axs[1].plot(averages[:,0], averages[:,6] - variances[:,6], linestyle='--', color='orange', linewidth=2)
# axs[1].plot(averages[:,0], averages[:,8], label="Dend_1-2", color='green', alpha=.5)
# axs[1].plot(averages[:,0], averages[:,8] + variances[:,8], linestyle='--', color='green', linewidth=2)
# axs[1].plot(averages[:,0], averages[:,8] - variances[:,8], linestyle='--', color='green', linewidth=2)
axs[1].set_ylabel(r'mRNAs', fontsize=20)



# axs[1].plot(num[:,0], num[:,3], color='red', zorder=1e4, label='Soma_theor')
# axs[1].plot(num[:,0], num[:,3] + num[:,4], linestyle='--', color='red')
# axs[1].plot(num[:,0], num[:,3] - num[:,4], linestyle='--', color='red')
# axs[1].plot(num[:,0], num[:,5], color='blue', zorder=1e4, label='Dend_1')
# axs[1].plot(num[:,0], num[:,5] + num[:,6], linestyle='--', color='blue')
# axs[1].plot(num[:,0], num[:,5] - num[:,6], linestyle='--', color='blue')
# axs[1].plot(num[:,0], num[:,7], color='orange', zorder=1e4, label='Dend_1-1')
# axs[1].plot(num[:,0], num[:,7] + num[:,8], linestyle='--', color='orange')
# axs[1].plot(num[:,0], num[:,7] - num[:,8], linestyle='--', color='orange')
# axs[1].plot(num[:,0], num[:,9], color='green', zorder=1e4, label='Dend_1-2')
# axs[1].plot(num[:,0], num[:,9] + num[:,10], linestyle='--', color='green')
# axs[1].plot(num[:,0], num[:,9] - num[:,10], linestyle='--', color='green')




############# Bulk proteins #############
# axs[2].plot(num[:,0], num[:,11], color='cyan', zorder=1e4, label='Soma_theor')
# axs[2].plot(num[:,0], num[:,11] + num[:,12], linestyle='--', color='cyan')
# axs[2].plot(num[:,0], num[:,11] - num[:,12], linestyle='--', color='cyan')
# axs[2].plot(num[:,0], num[:,13], color='pink', zorder=1e4, label='Dend_1_theor')
# axs[2].plot(num[:,0], num[:,13] + num[:,14], linestyle='--', color='pink')
# axs[2].plot(num[:,0], num[:,13] - num[:,14], linestyle='--', color='pink')
# axs[2].plot(num[:,0], num[:,19], color='brown', zorder=1e4, label='Dend_1_1_theor')
# axs[2].plot(num[:,0], num[:,19] + num[:,20], linestyle='--', color='brown')
# axs[2].plot(num[:,0], num[:,19] - num[:,20], linestyle='--', color='brown')

axs[2].plot(averages[:,0], averages[:,3], label="Soma", color='red', alpha=.7)
# axs[2].plot(num[:,0], averages[:,3] + variances[:,3], linestyle='--', color='red')
# axs[2].plot(num[:,0], averages[:,3] - variances[:,3], linestyle='--', color='red')

axs[2].plot(averages[:,0], averages[:,5], label="Dend_1", color='blue', alpha=.7)
axs[2].plot(averages[:,0], averages[:,5] + variances[:,5], linestyle='--', color='blue')
axs[2].plot(averages[:,0], averages[:,5] - variances[:,5], linestyle='--', color='blue')


# axs[2].plot(averages[:,0], averages[:,9], label="Dend_1-2", color='green', alpha=.7)
# axs[2].plot(averages[:,0], averages[:,9] + variances[:,9], linestyle='--', color='green')
# axs[2].plot(averages[:,0], averages[:,9] - variances[:,9], linestyle='--', color='green')

axs[2].set_ylabel(r'Bulk proteins', fontsize=20)



###################################################

axs[3].plot(num_exp[:,0], num_exp[:,6], label="Dend_1-1", color='red', alpha=.7)
axs[3].plot(averages[:,0], averages[:,6], label="Dend_1-1", color='blue', alpha=.7)
axs[3].plot(averages[:,0], averages[:,6] + variances[:,7], linestyle='--', color='blue')
axs[3].plot(averages[:,0], averages[:,6] - variances[:,7], linestyle='--', color='blue')


axs[3].plot(num_exp[:,0], num_exp[:,7], label="Dend_1-1", color='red', alpha=.7)
axs[3].plot(averages[:,0], averages[:,7], label="Dend_1-1", color='red', alpha=.7)
axs[3].plot(averages[:,0], averages[:,7] + variances[:,7], linestyle='--', color='red')
axs[3].plot(averages[:,0], averages[:,7] - variances[:,7], linestyle='--', color='red')




# axs[3].plot(num[:,0], num[:,15], color='maroon', zorder=1e4, label='Syn_1-1_theor')
# axs[3].plot(num[:,0], num[:,15] + num[:,16], linestyle='--', color='maroon')
# axs[3].plot(num[:,0], num[:,15] - num[:,16], linestyle='--', color='maroon')

# axs[3].plot(num[:,0], num[:,21], color='red', zorder=1e4, label='Syn_1-1_1_theor')
# axs[3].plot(num[:,0], num[:,21] + num[:,22], linestyle='--', color='red')
# axs[3].plot(num[:,0], num[:,21] - num[:,22], linestyle='--', color='red')

# axs[3].plot(num[:,0], num[:,29], color='red', zorder=1e4, label='Syn_1-2_2_theor')
# axs[3].plot(num[:,0], num[:,29] + num[:,30], linestyle='--', color='red')
# axs[3].plot(num[:,0], num[:,29] - num[:,30], linestyle='--', color='red')

axs[3].set_ylabel(r'Synaptic proteins', fontsize=18)



# axs[3].plot(averages[:,0], averages[:,10], label="Spine_1_1", color='maroon', alpha=.7)
# # axs[3].plot(num[:,0], averages[:,10] + variances[:,10], linestyle='--', color='maroon')
# # axs[3].plot(num[:,0], averages[:,10] - variances[:,10], linestyle='--', color='maroon')

# # axs[3].plot(averages[:,0], averages[:,12], label="Spine_1_2", alpha=.7)

# axs[3].plot(averages[:,0], averages[:,12], label="Spine_1-2_1", alpha=.5, color='red')
# # axs[3].plot(num[:,0], averages[:,12] + variances[:,12], linestyle='--', color='red')
# # axs[3].plot(num[:,0], averages[:,12] - variances[:,12], linestyle='--', color='red')

# axs[3].plot(averages[:,0], averages[:,14], label="Spine_1-2_2", alpha=.5, color='pink')
# # axs[3].plot(num[:,0], averages[:,14] + variances[:,14], linestyle='--', color='pink')
# # axs[3].plot(num[:,0], averages[:,14] - variances[:,14], linestyle='--', color='pink')

axs[3].set_xlabel(r'Time in hours', fontsize=18)

for ax in axs:
    ax.set_xlim([0,x_lim])
    ax.set_ylim(0)
    ax.legend(loc=4)

axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
