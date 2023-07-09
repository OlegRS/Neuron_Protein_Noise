import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

############ PARAMETERS #############
x_lim = 9999
sim_run_count = 1000 # Number of files with Gillespie simulations
mult_sim_run_count = sim_run_count
#####################################

fig_avrg, axs_avrg = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

# num = np.genfromtxt("../../data/nonstationary_expectations_and_std_deviations_100000_01_", delimiter=',')[:x_lim, :]
num = np.genfromtxt("../../data/soma_only/normal/normal", delimiter=',')[:x_lim, :]

# 1 - Active genes
# 2 - Soma mRNA
# 3 - Soma proteins - d_1 mRNA

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

# 1, soma__Gene
# 2
# 3, soma__mRNA
# 4
# 5, d_1__mRNA
# 6
# 7, d_1-1__mRNA
# 8
# 9, d_1-2__mRNA
# 10
# 11, soma__Prot
# 12
# 13, d_1__Prot
# 14
# 15, s_1_1__Prot
# 16
# 17, s_1_2__Prot
# 18
# 19, d_1-1__Prot
# 20
# 21, s_1-1_1__Prot
# 22
# 23, s_1-1_2__Prot
# 24
# 25, d_1-2__Prot
# 26
# 27, s_1-2_1__Prot
# 28
# 29, s_1-2_2__Prot
# 30

averages = np.zeros((x_lim, 4))#16))
variances = np.zeros((x_lim, 4))#16))

for file_id in range(sim_run_count):
    file_name = "../../data/soma_only/normal/sim/g_" + str(file_id)
    print(file_name)
    data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
    averages[:, 1:] += data[:, 1:]/mult_sim_run_count
    variances[:, 1:] += data[:, 1:]**2/mult_sim_run_count
#   axs_avrg[2].plot(data[:,0], data[:,15], label="Synapse_1-2_2", alpha=.5, color='pink')

# for file_id in range(sim_run_count):
#     file_name = "../../data/gillespie/test/1/g_" + str(file_id)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
#     averages[:, 1:] += data[:, 1:]/mult_sim_run_count

# for file_id in range(sim_run_count):
#     file_name = "../../data/gillespie/test/2/g_" + str(file_id)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
#     averages[:, 1:] += data[:, 1:]/mult_sim_run_count
#     variances[:, 1:] += data[:, 1:]**2/mult_sim_run_count

# for file_id in range(sim_run_count):
#     file_name = "../../data/gillespie/new/4/g_" + str(file_id)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
#     averages[:, 1:] += data[:, 1:]/mult_sim_run_count

averages[:, 0] = data[:, 0]
variances[:, 0] = data[:, 0]
variances[:, 1:] = np.sqrt(variances[:, 1:] - averages[:, 1:]**2)

# ################# Variances begin ####################
# for file_id in range(sim_run_count):
#     file_name = "../../data/gillespie/test/g_" + str(file_id)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
#     variances[:, 1:] += (data[:, 1:]-averages[:,1:])**2/mult_sim_run_count

# # for file_id in range(sim_run_count):
# #     file_name = "../../data/gillespie/new/2/g_" + str(file_id)
# #     print(file_name)
# #     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
# #     variances[:, 1:] += (data[:, 1:]-averages[:,1:])**2/mult_sim_run_count

# # for file_id in range(sim_run_count):
# #     file_name = "../../data/gillespie/new/3/g_" + str(file_id)
# #     print(file_name)
# #     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
# #     variances[:, 1:] += (data[:, 1:]-averages[:,1:])**2/mult_sim_run_count

# # for file_id in range(sim_run_count):
# #     file_name = "../../data/gillespie/new/4/g_" + str(file_id)
# #     print(file_name)
# #     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
# #     variances[:, 1:] += (data[:, 1:]-averages[:,1:])**2/mult_sim_run_count

# variances[:, 0] = data[:, 0]
# variances[:, 1:] = np.sqrt(variances[0:, 1:])

################# Variances end ####################

axs_avrg[0].plot(num[:,0], num[:,1], color='red', linewidth=2, zorder=1e4)
axs_avrg[0].plot(num[:,0], num[:,1] + num[:,2], linestyle='--', color='red', linewidth=2)
axs_avrg[0].plot(num[:,0], num[:,1] - num[:,2], linestyle='--', color='red', linewidth=2)


# axs_avrg[0].plot(num[:,0], num[:,3], color='red', linewidth=2, zorder=1e4)
# axs_avrg[0].plot(num[:,0], num[:,3] + num[:,4], linestyle='--', color='red', linewidth=2)
# axs_avrg[0].plot(num[:,0], num[:,3] - num[:,4], linestyle='--', color='red', linewidth=2)
# axs_avrg[0].plot(num[:,0], num[:,5], color='blue', linewidth=2, zorder=1e4)
# axs_avrg[0].plot(num[:,0], num[:,5] + num[:,6], linestyle='--', color='blue', linewidth=2)
# axs_avrg[0].plot(num[:,0], num[:,5] - num[:,6], linestyle='--', color='blue', linewidth=2)
# axs_avrg[0].plot(num[:,0], num[:,7], color='brown', linewidth=2, zorder=1e4)
# axs_avrg[0].plot(num[:,0], num[:,7] + num[:,8], linestyle='--', color='brown', linewidth=2)
# axs_avrg[0].plot(num[:,0], num[:,7] - num[:,8], linestyle='--', color='brown', linewidth=2)
axs_avrg[0].set_ylabel(r'mRNA count', fontsize=20)


axs_avrg[0].plot(averages[:,0], averages[:,1], label='Soma', color='red', alpha=.5)
axs_avrg[0].plot(averages[:,0], averages[:,1] + variances[:,1], linestyle='--', color='red', linewidth=2)
axs_avrg[0].plot(averages[:,0], averages[:,1] - variances[:,1], linestyle='--', color='red', linewidth=2)


# axs_avrg[0].plot(averages[:,0], averages[:,2], label='Soma', color='red', alpha=.5)
# axs_avrg[0].plot(averages[:,0], averages[:,2] + variances[:,2], linestyle='--', color='red', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,2] - variances[:,2], linestyle='--', color='red', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,4], label="Dend_1", color='blue', alpha=.5)
# axs_avrg[0].plot(averages[:,0], averages[:,4] + variances[:,4], linestyle='--', color='blue', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,4] - variances[:,4], linestyle='--', color='blue', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,6], label="Dend_1-1", color='orange', alpha=.5)
# axs_avrg[0].plot(averages[:,0], averages[:,6] + variances[:,6], linestyle='--', color='orange', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,6] - variances[:,6], linestyle='--', color='orange', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,8], label="Dend_1-2", color='green', alpha=.5)
# axs_avrg[0].plot(averages[:,0], averages[:,8] + variances[:,8], linestyle='--', color='green', linewidth=2)
# axs_avrg[0].plot(averages[:,0], averages[:,8] - variances[:,8], linestyle='--', color='green', linewidth=2)


axs_avrg[1].plot(num[:,0], num[:,3], color='cyan', zorder=1e4, label='Soma_theor')
axs_avrg[1].plot(num[:,0], num[:,3] + num[:,4], linestyle='--', color='cyan')
axs_avrg[1].plot(num[:,0], num[:,3] - num[:,4], linestyle='--', color='cyan')

# axs_avrg[1].plot(num[:,0], num[:,11], color='cyan', zorder=1e4, label='Soma_theor')
# axs_avrg[1].plot(num[:,0], num[:,11] + num[:,12], linestyle='--', color='cyan')
# axs_avrg[1].plot(num[:,0], num[:,11] - num[:,12], linestyle='--', color='cyan')
# axs_avrg[1].plot(num[:,0], num[:,13], color='pink', zorder=1e4, label='Dend_1_theor')
# axs_avrg[1].plot(num[:,0], num[:,13] + num[:,14], linestyle='--', color='pink')
# axs_avrg[1].plot(num[:,0], num[:,13] - num[:,14], linestyle='--', color='pink')
# axs_avrg[1].plot(num[:,0], num[:,19], color='brown', zorder=1e4, label='Dend_1_1_theor')
# axs_avrg[1].plot(num[:,0], num[:,19] + num[:,20], linestyle='--', color='brown')
# axs_avrg[1].plot(num[:,0], num[:,19] - num[:,20], linestyle='--', color='brown')
axs_avrg[1].set_ylabel(r'Dendritic protein count', fontsize=20)

axs_avrg[1].plot(averages[:,0], averages[:,2], label="Soma mRNA", color='blue', alpha=.7)
axs_avrg[1].plot(num[:,0], averages[:,2] + variances[:,2], linestyle='--', color='cyan')
axs_avrg[1].plot(num[:,0], averages[:,2] - variances[:,2], linestyle='--', color='cyan')


# axs_avrg[1].plot(averages[:,0], averages[:,3], label="Soma", color='blue', alpha=.7)
# axs_avrg[1].plot(num[:,0], averages[:,3] + variances[:,3], linestyle='--', color='cyan')
# axs_avrg[1].plot(num[:,0], averages[:,3] - variances[:,3], linestyle='--', color='cyan')
# axs_avrg[1].plot(averages[:,0], averages[:,5], label="Dend_1", color='red', alpha=.7)
# axs_avrg[1].plot(averages[:,0], averages[:,7], label="Dend_1-1", color='orange', alpha=.7)
# axs_avrg[1].plot(averages[:,0], averages[:,9], label="Dend_1-2", color='green', alpha=.7)
# axs_avrg[1].plot(averages[:,0], averages[:,9] + variances[:,9], linestyle='--', color='green')
# axs_avrg[1].plot(averages[:,0], averages[:,9] - variances[:,9], linestyle='--', color='green')


axs_avrg[2].plot(num[:,0], num[:,5], color='red', zorder=1e4, label='Syn_1-1_theor')
axs_avrg[2].plot(num[:,0], num[:,5] + num[:,6], linestyle='--', color='red')
axs_avrg[2].plot(num[:,0], num[:,5] - num[:,6], linestyle='--', color='red')


# axs_avrg[2].plot(num[:,0], num[:,15], color='maroon', zorder=1e4, label='Syn_1-1_theor')
# axs_avrg[2].plot(num[:,0], num[:,15] + num[:,16], linestyle='--', color='maroon')
# axs_avrg[2].plot(num[:,0], num[:,15] - num[:,16], linestyle='--', color='maroon')
# axs_avrg[2].plot(num[:,0], num[:,21], color='red', zorder=1e4, label='Syn_1-1_1_theor')
# axs_avrg[2].plot(num[:,0], num[:,21] + num[:,22], linestyle='--', color='red')
# axs_avrg[2].plot(num[:,0], num[:,21] - num[:,22], linestyle='--', color='red')
# axs_avrg[2].plot(num[:,0], num[:,29], color='red', zorder=1e4, label='Syn_1-1_1_theor')
# axs_avrg[2].plot(num[:,0], num[:,29] + num[:,30], linestyle='--', color='red')
# axs_avrg[2].plot(num[:,0], num[:,29] - num[:,30], linestyle='--', color='red')

axs_avrg[2].set_ylabel(r'Synaptic protein count', fontsize=20)
axs_avrg[2].set_xlabel(r'Time in hours', fontsize=20)


axs_avrg[2].plot(averages[:,0], averages[:,3], label="Soma_1_1", color='maroon', alpha=.7)
axs_avrg[2].plot(num[:,0], averages[:,3] + variances[:,3], linestyle='--', color='maroon')
axs_avrg[2].plot(num[:,0], averages[:,3] - variances[:,3], linestyle='--', color='maroon')


# axs_avrg[2].plot(averages[:,0], averages[:,10], label="Synapse_1_1", color='maroon', alpha=.7)
# axs_avrg[2].plot(num[:,0], averages[:,10] + variances[:,10], linestyle='--', color='maroon')
# axs_avrg[2].plot(num[:,0], averages[:,10] - variances[:,10], linestyle='--', color='maroon')
# axs_avrg[2].plot(averages[:,0], averages[:,11], label="Synapse_1_2", alpha=.7)
# axs_avrg[2].plot(averages[:,0], averages[:,12], label="Synapse_1-1_1", alpha=.7)
# axs_avrg[2].plot(averages[:,0], averages[:,13], label="Synapse_1-1_2", alpha=.7)
# axs_avrg[2].plot(averages[:,0], averages[:,14], label="Synapse_1-2_1", alpha=.5, color='red')
# axs_avrg[2].plot(num[:,0], averages[:,14] + variances[:,14], linestyle='--', color='red')
# axs_avrg[2].plot(num[:,0], averages[:,14] - variances[:,14], linestyle='--', color='red')
# axs_avrg[2].plot(averages[:,0], averages[:,15], label="Synapse_1-2_2", alpha=.5, color='pink')
# axs_avrg[2].plot(num[:,0], averages[:,15] + variances[:,15], linestyle='--', color='pink')
# axs_avrg[2].plot(num[:,0], averages[:,15] - variances[:,15], linestyle='--', color='pink')

for ax in axs_avrg:
    ax.set_xlim([0,x_lim])
    ax.set_ylim(0)
    ax.legend(loc=4)

# axs_avrg[2].set_ylim([0,10000])

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)

fig_var, axs_var = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))


axs_var[0].plot(averages[:,0], variances[:,1], linestyle='--', color='green')
axs_var[0].plot(num[:,0], num[:, 2], linestyle='-', color='red')
axs_var[0].set_ylabel(r'Active genes $\sigma$', fontsize=20)

axs_var[1].plot(variances[:,0], variances[:,2], linestyle='--', color='pink')
axs_var[1].plot(num[:,0], num[:,4], linestyle='--', color='red')
axs_var[1].set_ylabel(r'mRNAs $\sigma$', fontsize=20)

axs_var[2].plot(variances[:,0], variances[:,3], linestyle='--', color='pink')
axs_var[2].plot(num[:,0], num[:,6], linestyle='--', color='red')
axs_var[2].set_ylabel(r'Proteins $\sigma$', fontsize=20)
axs_var[2].set_xlabel('Time in hours', fontsize=20)


axs_var[0].plot(results[0]['time'], np.sqrt(gene_var))
axs_var[0].grid()
axs_var[1].plot(results[0]['time'], np.sqrt(mRNA_var))
axs_var[1].grid()
axs_var[2].plot(results[0]['time'], np.sqrt(prot_var))
axs_var[2].grid()


axs_avrg[0].plot(results[0]['time'], gene_average)
axs_avrg[0].grid()
axs_avrg[1].plot(results[0]['time'], mRNA_average)
axs_avrg[1].grid()
axs_avrg[2].plot(results[0]['time'], prot_average)
axs_avrg[2].grid()


plt.show()
