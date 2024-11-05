# averages_1 = np.zeros((n_points, n_compartments))
# variances_1 = np.zeros((n_points, n_compartments))
# for file_id in range(sim_run_count_1):
#     print("Loading file " + str(file_id+2))
#     file_name = "../../../data/gillespie/stationary_moments_for_cosyne_2_genes_more_more/SM_" + str(file_id+2)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
#     averages_1[:, 1:] += data[:, 1:]/sim_run_count_1
#     variances_1[:, 1:] += data[:, 1:]**2/sim_run_count_1


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(18/2*1.2, 3.2*1.9/2*1.2))

####################### COL_1 ###########################
axs[0].axhline(num_exp[5], color='#4CFFFF', linewidth=3, zorder=10, label='S_1_1') # s_1_1 prot
axs[0].axhline(num_exp[5] + num_std[5], linestyle='--', color='#4CFFFF', linewidth=3, zorder=10, label='S_1_1 $\pm\sigma$')
axs[0].axhline(num_exp[5] - num_std[5], linestyle='--', color='#4CFFFF', linewidth=3, zorder=10)
axs[0].plot(data[:,0], averages_1[:,10], color='#4CFFFF', alpha=.8)

axs[0].axhline(num_exp[6], color='#4C4CA6', linewidth=3, zorder=10, label='S_1_2') # s_1_2 prot
axs[0].axhline(num_exp[6] + num_std[6], linestyle='--', color='#4C4CA6', linewidth=3, zorder=10, label='S_1_2 $\pm\sigma$')
axs[0].axhline(num_exp[6] - num_std[6], linestyle='--', color='#4C4CA6', linewidth=3, zorder=10)
axs[0].plot(data[:,0], averages_1[:,11], color='#4C4CA6', alpha=.8)

axs[0].axhline(num_exp[9], color='#FA8072', linewidth=3, zorder=10, label='D_1-1_1') # s_1-1_1 prot
axs[0].axhline(num_exp[9] + num_std[9], linestyle='--', color='#FA8072', linewidth=3, zorder=10, label='D_1-1_1 $\pm\sigma$')
axs[0].axhline(num_exp[9] - num_std[9], linestyle='--', color='#FA8072', linewidth=3, zorder=10)
axs[0].plot(data[:,0], averages_1[:,12], alpha=.8, color='#FA8072')

axs[0].axhline(num_exp[10], color='#FFD700', linewidth=3, zorder=10, label='S_1-1_2') # s_1-1_2 prot
axs[0].axhline(num_exp[10] + num_std[10], linestyle='--', color='#FFD700', linewidth=3, zorder=10, label='S_1-1_2 $\pm\sigma$')
axs[0].axhline(num_exp[10] - num_std[10], linestyle='--', color='#FFD700', linewidth=3, zorder=10)
axs[0].plot(data[:,0], averages_1[:,13], alpha=.8, color='#FFD700')

axs[0].axhline(num_exp[13], color='#7CFC00', linewidth=3, zorder=10, label='S_1-2_1') # s_1-2_1 prot
axs[0].axhline(num_exp[13] + num_std[13], linestyle='--', color='#7CFC00', linewidth=3, zorder=10, label='S_1-2_1 $\pm\sigma$')
axs[0].axhline(num_exp[13] - num_std[13], linestyle='--', color='#7CFC00', linewidth=3, zorder=10)
axs[0].plot(data[:,0], averages_1[:,14], alpha=.8, color='#7CFC00')

axs[0].axhline(num_exp[14], color='#228B22', linewidth=3, zorder=10, label='S_1-2_2') # s_1-2_2 prot
axs[0].axhline(num_exp[14] + num_std[14], linestyle='--', color='#228B22', linewidth=3, zorder=10, label='S_1-2_2 $\pm\sigma$')
axs[0].axhline(num_exp[14] - num_std[14], linestyle='--', color='#228B22', linewidth=3, zorder=10)
axs[0].plot(data[:,0], averages_1[:,15], alpha=.8, color='#228B22')

axs[0].set_facecolor('white')

####################### COL_2 ###########################
axs[1].axhline(num_exp[5], color='#4CFFFF', linewidth=3) # s_1_1 prot
axs[1].axhline(num_exp[5] + num_std[5], linestyle='--', color='#4CFFFF', linewidth=3)
axs[1].axhline(num_exp[5] - num_std[5], linestyle='--', color='#4CFFFF', linewidth=3)
axs[1].plot(data[:,0], averages_2[:,10], label='S_1_1', color='#4CFFFF', alpha=.8)
axs[1].plot(data[:,0], averages_2[:,10] + stds_2[:,10], label='S_1_1 $\pm$ SD', alpha=.5, color='#4CFFFF')
axs[1].plot(data[:,0], averages_2[:,10] - stds_2[:,10], alpha=.5, color='#4CFFFF')

axs[1].axhline(num_exp[6], color='#4C4CA6', linewidth=3) # s_1_2 prot
axs[1].axhline(num_exp[6] + num_std[6], linestyle='--', color='#4C4CA6', linewidth=3)
axs[1].axhline(num_exp[6] - num_std[6], linestyle='--', color='#4C4CA6', linewidth=3)
axs[1].plot(data[:,0], averages_2[:,11], label='S_1_2', alpha=.8, color='#4C4CA6')
axs[1].plot(data[:,0], averages_2[:,11] + stds_2[:,11], label='S_1_2 $\pm$ SD', alpha=.5, color='#4C4CA6')
axs[1].plot(data[:,0], averages_2[:,11] - stds_2[:,11], alpha=.5, color='#4C4CA6')

axs[1].axhline(num_exp[9], color='#FA8072', linewidth=3) # s_1-1_1 prot
axs[1].axhline(num_exp[9] + num_std[9], linestyle='--', color='#FA8072', linewidth=3)
axs[1].axhline(num_exp[9] - num_std[9], linestyle='--', color='#FA8072', linewidth=3)
axs[1].plot(data[:,0], averages_2[:,12], label='S_1-1_1', alpha=.8, color='#FA8072')
axs[1].plot(data[:,0], averages_2[:,12] + stds_2[:,12], label='S_1-1_1 $\pm$ SD', alpha=.5, color='#FA8072')
axs[1].plot(data[:,0], averages_2[:,12] - stds_2[:,12], alpha=.5, color='#FA8072')

axs[1].axhline(num_exp[10], color='#FFD700', linewidth=3) # s_1-1_2 prot
axs[1].axhline(num_exp[10] + num_std[10], linestyle='--', color='#FFD700', linewidth=3)
axs[1].axhline(num_exp[10] - num_std[10], linestyle='--', color='#FFD700', linewidth=3)
axs[1].plot(data[:,0], averages_2[:,13], label='S_1-1_2', alpha=.8, color='#FFD700')
axs[1].plot(data[:,0], averages_2[:,13] + stds_2[:,13], label='S_1-1_2 $\pm$ SD', alpha=.5, color='#FFD700')
axs[1].plot(data[:,0], averages_2[:,13] - stds_2[:,13], alpha=.5, color='#FFD700')

axs[1].axhline(num_exp[13], color='#7CFC00', linewidth=3) # s_1-2_1 prot
axs[1].axhline(num_exp[13] + num_std[13], linestyle='--', color='#7CFC00', linewidth=3)
axs[1].axhline(num_exp[13] - num_std[13], linestyle='--', color='#7CFC00', linewidth=3)
axs[1].plot(data[:,0], averages_2[:,14], label='S_1-2_1', alpha=.8, color='#7CFC00')
axs[1].plot(data[:,0], averages_2[:,14] + stds_2[:,14], label='S_1-2_1 $\pm$ SD', alpha=.5, color='#7CFC00')
axs[1].plot(data[:,0], averages_2[:,14] - stds_2[:,14], alpha=.5, color='#7CFC00')

axs[1].axhline(num_exp[14], color='#228B22', linewidth=3) # s_1-2_2 prot
axs[1].axhline(num_exp[14] + num_std[14], linestyle='--', color='#228B22', linewidth=3)
axs[1].axhline(num_exp[14] - num_std[14], linestyle='--', color='#228B22', linewidth=3)
axs[1].plot(data[:,0], averages_2[:,15], label='S_1-2_2', alpha=.8, color='#228B22')
axs[1].plot(data[:,0], averages_2[:,15] + stds_2[:,15], label='S_1-2_2 $\pm$ SD', alpha=.5, color='#228B22')
axs[1].plot(data[:,0], averages_2[:,15] - stds_2[:,15], alpha=.5, color='#228B22')

axs[1].set_facecolor('white')

for i in range(2):
    axs[i].set_xlabel(r'Time in hours', fontsize=12)
    # ax.set_xlim([0,num_exp_1[num_exp_1.shape[0]-1,0]])
    axs[i].set_xlim([0,1000])
    axs[i].set_ylim(0)
    # axs[i].legend(loc=4, fontsize=6.15, handlelength=3).set_zorder(100)
    axs[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# axs[0].set_title("Average counts over multiple trajectories", fontsize=12)

axs[0].set_ylabel('Synaptic protein count', fontsize=12)


axs[0].set_ylim([0, 5.6e3])
axs[1].set_ylim([0, 5.6e3])




plt.tight_layout()

fig.patch.set_facecolor('none')
plt.savefig('/home/oleg/oleg_windows/olegr/Desktop/study/BBSRC_application/stationary_num_MC_2_genes_new.png', dpi=300, bbox_inches='tight')
# plt.savefig('/home/oleg/sync/study/ulster/BBSRC_representational_drift/img/stationary_num_MC_2_genes_new.pdf', dpi=300)
# plt.savefig('/home/oleg/sync/study/ulster/BBSRC_representational_drift/img/stationary_num_MC_2_genes_new.png', dpi=300)

plt.show()

fig.clear()
