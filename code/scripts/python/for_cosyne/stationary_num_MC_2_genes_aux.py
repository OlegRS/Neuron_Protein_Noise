# averages_1 = np.zeros((n_points, n_compartments))
# variances_1 = np.zeros((n_points, n_compartments))
# for file_id in range(sim_run_count_1):
#     print("Loading file " + str(file_id+2))
#     file_name = "../../../data/gillespie/stationary_moments_for_cosyne_2_genes_more_more/SM_" + str(file_id+2)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[:, 1:]
#     averages_1[:, 1:] += data[:, 1:]/sim_run_count_1
#     variances_1[:, 1:] += data[:, 1:]**2/sim_run_count_1


fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(9, 3.2*1.9*1.8))

####################### COL_1 ###########################
axs[0,0].set_title('Single simulation', fontsize=15)
axs[0,1].set_title('Averaged 500 simulations', fontsize=15)

axs[0,0].axhline(num_exp[0], color='#FF4C4C', linewidth=3, label='Soma', zorder=10) # Soma mRNAs
axs[0,0].axhline(num_exp[0] + num_std[0], linestyle='--', color='#FF4C4C', linewidth=3, label='Soma $\pm\sigma$', zorder=10)
axs[0,0].axhline(num_exp[0] - num_std[0], linestyle='--', color='#FF4C4C', linewidth=3, zorder=10)
axs[0,0].plot(data[:,0], averages_1[:,1], color='#FF4C4C', alpha=.8)
axs[0,0].set_facecolor('white')


axs[1,0].axhline(num_exp[1], color='#FF4C4C', linewidth=3, label='Soma', zorder=10) # Soma mRNAs
axs[1,0].axhline(num_exp[1] + num_std[1], linestyle='--', color='#FF4C4C', linewidth=3, label='Soma $\pm\sigma$', zorder=10)
axs[1,0].axhline(num_exp[1] - num_std[1], linestyle='--', color='#FF4C4C', linewidth=3, zorder=10)
axs[1,0].plot(data[:,0], averages_1[:,2], color='#FF4C4C', alpha=.8)

axs[1,0].axhline(num_exp[3], color='#4C4CFF', linewidth=3, label='D_1', zorder=10) # d_1 mRNAs
axs[1,0].axhline(num_exp[3] + num_std[3], linestyle='--', color='#4C4CFF', linewidth=3, label='D_1 $\pm\sigma$', zorder=10)
axs[1,0].axhline(num_exp[3] - num_std[3], linestyle='--', color='#4C4CFF', linewidth=3, zorder=10)
axs[1,0].plot(data[:,0], averages_1[:,4], alpha=.8, color='#4C4CFF')

axs[1,0].axhline(num_exp[7], color='#FFBF4C', linewidth=3, label='D_1-1', zorder=10) # d_1_1 mRNAs
axs[1,0].axhline(num_exp[7] + num_std[7], linestyle='--', color='#FFBF4C', linewidth=3, label='D_1-1 $\pm\sigma$', zorder=10)
axs[1,0].axhline(num_exp[7] - num_std[7], linestyle='--', color='#FFBF4C', linewidth=3, zorder=10)
axs[1,0].plot(data[:,0], averages_1[:,6], alpha=.8, color='#FFBF4C')

axs[1,0].axhline(num_exp[11], color='#4CA64C', linewidth=3, label='D_1-2', zorder=10) # d_1_2 mRNAs
axs[1,0].axhline(num_exp[11] + num_std[11], linestyle='--', color='#4CA64C', linewidth=3, label='D_1-2 $\pm\sigma$', zorder=10)
axs[1,0].axhline(num_exp[11] - num_std[11], linestyle='--', color='#4CA64C', linewidth=3, zorder=10)
axs[1,0].plot(data[:,0], averages_1[:,8], alpha=.8, color='#4CA64C')

axs[1,0].set_facecolor('white')


axs[2,0].axhline(num_exp[2], color='#FF4C4C', linewidth=3, label='Soma', zorder=10) # Soma prot
axs[2,0].axhline(num_exp[2] + num_std[2], linestyle='--', color='#FF4C4C', linewidth=3, label='Soma $\pm\sigma$', zorder=10)
axs[2,0].axhline(num_exp[2] - num_std[2], linestyle='--', color='#FF4C4C', linewidth=3, zorder=10)
axs[2,0].plot(data[:,0], averages_1[:,3], color='#FF4C4C', alpha=.8)

axs[2,0].axhline(num_exp[4], color='#4C4CFF', linewidth=3, zorder=10, label='D_1') # d_1 prot
axs[2,0].axhline(num_exp[4] + num_std[4], linestyle='--', color='#4C4CFF', linewidth=3, zorder=10, label='D_1 $\pm\sigma$')
axs[2,0].axhline(num_exp[4] - num_std[4], linestyle='--', color='#4C4CFF', linewidth=3, zorder=10)
axs[2,0].plot(data[:,0], averages_1[:,5], color='#4C4CFF', alpha=.8)

axs[2,0].axhline(num_exp[8], color='#FFBF4C', linewidth=3, zorder=10, label='D_1-1') # d_1_1 prot
axs[2,0].axhline(num_exp[8] + num_std[8], linestyle='--', color='#FFBF4C', linewidth=3, zorder=10, label='D_1-1 $\pm\sigma$')
axs[2,0].axhline(num_exp[8] - num_std[8], linestyle='--', color='#FFBF4C', linewidth=3, zorder=10)
axs[2,0].plot(data[:,0], averages_1[:,7], alpha=.8, color='#FFBF4C')

axs[2,0].axhline(num_exp[12], color='#4CA64C', linewidth=3, zorder=10, label='D_1-2') # d_1_2 prot
axs[2,0].axhline(num_exp[12] + num_std[12], linestyle='--', color='#4CA64C', linewidth=3, zorder=10, label='D_1-2 $\pm\sigma$')
axs[2,0].axhline(num_exp[12] - num_std[12], linestyle='--', color='#4CA64C', linewidth=3, zorder=10)
axs[2,0].plot(data[:,0], averages_1[:,9], alpha=.8, color='#4CA64C')

axs[2,0].set_facecolor('white')


axs[3,0].axhline(num_exp[5], color='#4CFFFF', linewidth=3, zorder=10, label='S_1_1') # s_1_1 prot
axs[3,0].axhline(num_exp[5] + num_std[5], linestyle='--', color='#4CFFFF', linewidth=3, zorder=10, label='S_1_1 $\pm\sigma$')
axs[3,0].axhline(num_exp[5] - num_std[5], linestyle='--', color='#4CFFFF', linewidth=3, zorder=10)
axs[3,0].plot(data[:,0], averages_1[:,10], color='#4CFFFF', alpha=.8)

axs[3,0].axhline(num_exp[6], color='#4C4CA6', linewidth=3, zorder=10, label='S_1_2') # s_1_2 prot
axs[3,0].axhline(num_exp[6] + num_std[6], linestyle='--', color='#4C4CA6', linewidth=3, zorder=10, label='S_1_2 $\pm\sigma$')
axs[3,0].axhline(num_exp[6] - num_std[6], linestyle='--', color='#4C4CA6', linewidth=3, zorder=10)
axs[3,0].plot(data[:,0], averages_1[:,11], color='#4C4CA6', alpha=.8)

axs[3,0].axhline(num_exp[9], color='#FA8072', linewidth=3, zorder=10, label='D_1-1_1') # s_1-1_1 prot
axs[3,0].axhline(num_exp[9] + num_std[9], linestyle='--', color='#FA8072', linewidth=3, zorder=10, label='D_1-1_1 $\pm\sigma$')
axs[3,0].axhline(num_exp[9] - num_std[9], linestyle='--', color='#FA8072', linewidth=3, zorder=10)
axs[3,0].plot(data[:,0], averages_1[:,12], alpha=.8, color='#FA8072')

axs[3,0].axhline(num_exp[10], color='#FFD700', linewidth=3, zorder=10, label='S_1-1_2') # s_1-1_2 prot
axs[3,0].axhline(num_exp[10] + num_std[10], linestyle='--', color='#FFD700', linewidth=3, zorder=10, label='S_1-1_2 $\pm\sigma$')
axs[3,0].axhline(num_exp[10] - num_std[10], linestyle='--', color='#FFD700', linewidth=3, zorder=10)
axs[3,0].plot(data[:,0], averages_1[:,13], alpha=.8, color='#FFD700')

axs[3,0].axhline(num_exp[13], color='#7CFC00', linewidth=3, zorder=10, label='S_1-2_1') # s_1-2_1 prot
axs[3,0].axhline(num_exp[13] + num_std[13], linestyle='--', color='#7CFC00', linewidth=3, zorder=10, label='S_1-2_1 $\pm\sigma$')
axs[3,0].axhline(num_exp[13] - num_std[13], linestyle='--', color='#7CFC00', linewidth=3, zorder=10)
axs[3,0].plot(data[:,0], averages_1[:,14], alpha=.8, color='#7CFC00')

axs[3,0].axhline(num_exp[14], color='#228B22', linewidth=3, zorder=10, label='S_1-2_2') # s_1-2_2 prot
axs[3,0].axhline(num_exp[14] + num_std[14], linestyle='--', color='#228B22', linewidth=3, zorder=10, label='S_1-2_2 $\pm\sigma$')
axs[3,0].axhline(num_exp[14] - num_std[14], linestyle='--', color='#228B22', linewidth=3, zorder=10)
axs[3,0].plot(data[:,0], averages_1[:,15], alpha=.8, color='#228B22')

axs[3,0].set_facecolor('white')

####################### COL_2 ###########################
axs[0,1].axhline(num_exp[0], color='#FF4C4C', linewidth=3) # Soma mRNAs
axs[0,1].axhline(num_exp[0] + num_std[0], linestyle='--', color='#FF4C4C', linewidth=3)
axs[0,1].axhline(num_exp[0] - num_std[0], linestyle='--', color='#FF4C4C', linewidth=3)
axs[0,1].plot(data[:,0], averages_2[:,1], label='Soma', color='#FF4C4C', alpha=.8)
axs[0,1].plot(data[:,0], averages_2[:,1] + stds_2[:,1], label='Soma $\pm$ SD', color='#FF4C4C', alpha=.3)
axs[0,1].plot(data[:,0], averages_2[:,1] - stds_2[:,1], color='#FF4C4C', alpha=.3)

axs[0,1].set_facecolor('white')


axs[1,1].axhline(num_exp[1], color='#FF4C4C', linewidth=3) # Soma mRNAs
axs[1,1].axhline(num_exp[1] + num_std[1], linestyle='--', color='#FF4C4C', linewidth=3)
axs[1,1].axhline(num_exp[1] - num_std[1], linestyle='--', color='#FF4C4C', linewidth=3)
axs[1,1].plot(data[:,0], averages_2[:,2], label='Soma', color='#FF4C4C', alpha=.8)
axs[1,1].plot(data[:,0], averages_2[:,2] + stds_2[:,2], label='Soma $\pm$ SD', color='#FF4C4C', alpha=.3)
axs[1,1].plot(data[:,0], averages_2[:,2] - stds_2[:,2], color='#FF4C4C', alpha=.3)

axs[1,1].axhline(num_exp[3], color='#4C4CFF', linewidth=3) # d_1 mRNAs
axs[1,1].axhline(num_exp[3] + num_std[3], linestyle='--', color='#4C4CFF', linewidth=3)
axs[1,1].axhline(num_exp[3] - num_std[3], linestyle='--', color='#4C4CFF', linewidth=3)
axs[1,1].plot(data[:,0], averages_2[:,4], label='D_1', alpha=.8, color='#4C4CFF')
axs[1,1].plot(data[:,0], averages_2[:,4] + stds_2[:,4], label='D_1 $\pm$ SD', alpha=.3, color='#4C4CFF')
axs[1,1].plot(data[:,0], averages_2[:,4] - stds_2[:,4], alpha=.3, color='#4C4CFF')

axs[1,1].axhline(num_exp[7], color='#FFBF4C', linewidth=3) # d_1_1 mRNAs
axs[1,1].axhline(num_exp[7] + num_std[7], linestyle='--', color='#FFBF4C', linewidth=3)
axs[1,1].axhline(num_exp[7] - num_std[7], linestyle='--', color='#FFBF4C', linewidth=3)
axs[1,1].plot(data[:,0], averages_2[:,6], label='D_1-1', alpha=.8, color='#FFBF4C')
axs[1,1].plot(data[:,0], averages_2[:,6] + stds_2[:,6], label='D_1-1 $\pm$ SD', alpha=.3, color='#FFBF4C')
axs[1,1].plot(data[:,0], averages_2[:,6] - stds_2[:,6], alpha=.3, color='#FFBF4C')

axs[1,1].axhline(num_exp[11], color='#4CA64C', linewidth=3) # d_1_2 mRNAs
axs[1,1].axhline(num_exp[11] + num_std[11], linestyle='--', color='#4CA64C', linewidth=3)
axs[1,1].axhline(num_exp[11] - num_std[11], linestyle='--', color='#4CA64C', linewidth=3)
axs[1,1].plot(data[:,0], averages_2[:,8], label='D_1-2', alpha=.8, color='#4CA64C')
axs[1,1].plot(data[:,0], averages_2[:,8] + stds_2[:,8], label='D_1-2 $\pm$ SD', alpha=.3, color='#4CA64C')
axs[1,1].plot(data[:,0], averages_2[:,8] - stds_2[:,8], alpha=.3, color='#4CA64C')

axs[1,1].set_facecolor('white')


axs[2,1].axhline(num_exp[2], color='#FF4C4C', linewidth=3) # Soma prot
axs[2,1].axhline(num_exp[2] + num_std[2], linestyle='--', color='#FF4C4C', linewidth=3)
axs[2,1].axhline(num_exp[2] - num_std[2], linestyle='--', color='#FF4C4C', linewidth=3)
axs[2,1].plot(data[:,0], averages_2[:,3], label='Soma', color='#FF4C4C', alpha=.8)
axs[2,1].plot(data[:,0], averages_2[:,3] + stds_2[:,3], label='Soma $\pm$ SD', alpha=.3, color='#FF4C4C')
axs[2,1].plot(data[:,0], averages_2[:,3] - stds_2[:,3], alpha=.3, color='#FF4C4C')

axs[2,1].axhline(num_exp[4], color='#4C4CFF', linewidth=3) # d_1 prot
axs[2,1].axhline(num_exp[4] + num_std[4], linestyle='--', color='#4C4CFF', linewidth=3)
axs[2,1].axhline(num_exp[4] - num_std[4], linestyle='--', color='#4C4CFF', linewidth=3)
axs[2,1].plot(data[:,0], averages_2[:,5], label='D_1', alpha=.8, color='#4C4CFF')
axs[2,1].plot(data[:,0], averages_2[:,5] + stds_2[:,5], label='D_1 $\pm$ SD', alpha=.3, color='#4C4CFF')
axs[2,1].plot(data[:,0], averages_2[:,5] - stds_2[:,5], alpha=.3, color='#4C4CFF')

axs[2,1].axhline(num_exp[8], color='#FFBF4C', linewidth=3) # d_1_1 prot
axs[2,1].axhline(num_exp[8] + num_std[8], linestyle='--', color='#FFBF4C', linewidth=3)
axs[2,1].axhline(num_exp[8] - num_std[8], linestyle='--', color='#FFBF4C', linewidth=3)
axs[2,1].plot(data[:,0], averages_2[:,7], label='D_1-1', alpha=.8, color='#FFBF4C')
axs[2,1].plot(data[:,0], averages_2[:,7] + stds_2[:,7], label='D_1-1 $\pm$ SD', alpha=.3, color='#FFBF4C')
axs[2,1].plot(data[:,0], averages_2[:,7] - stds_2[:,7], alpha=.3, color='#FFBF4C')

axs[2,1].axhline(num_exp[12], color='#4CA64C', linewidth=3) # d_1_2 prot
axs[2,1].axhline(num_exp[12] + num_std[12], linestyle='--', color='#4CA64C', linewidth=3)
axs[2,1].axhline(num_exp[12] - num_std[12], linestyle='--', color='#4CA64C', linewidth=3)
axs[2,1].plot(data[:,0], averages_2[:,9], label='D_1-2', alpha=.8, color='#4CA64C')
axs[2,1].plot(data[:,0], averages_2[:,9] + stds_2[:,9], label='D_1-2 $\pm$ SD', alpha=.3, color='#4CA64C')
axs[2,1].plot(data[:,0], averages_2[:,9] - stds_2[:,9], alpha=.3, color='#4CA64C')

axs[2,1].set_facecolor('white')


axs[3,1].axhline(num_exp[5], color='#4CFFFF', linewidth=3) # s_1_1 prot
axs[3,1].axhline(num_exp[5] + num_std[5], linestyle='--', color='#4CFFFF', linewidth=3)
axs[3,1].axhline(num_exp[5] - num_std[5], linestyle='--', color='#4CFFFF', linewidth=3)
axs[3,1].plot(data[:,0], averages_2[:,10], label='S_1_1', color='#4CFFFF', alpha=.8)
axs[3,1].plot(data[:,0], averages_2[:,10] + stds_2[:,10], label='S_1_1 $\pm$ SD', alpha=.5, color='#4CFFFF')
axs[3,1].plot(data[:,0], averages_2[:,10] - stds_2[:,10], alpha=.5, color='#4CFFFF')

axs[3,1].axhline(num_exp[6], color='#4C4CA6', linewidth=3) # s_1_2 prot
axs[3,1].axhline(num_exp[6] + num_std[6], linestyle='--', color='#4C4CA6', linewidth=3)
axs[3,1].axhline(num_exp[6] - num_std[6], linestyle='--', color='#4C4CA6', linewidth=3)
axs[3,1].plot(data[:,0], averages_2[:,11], label='S_1_2', alpha=.8, color='#4C4CA6')
axs[3,1].plot(data[:,0], averages_2[:,11] + stds_2[:,11], label='S_1_2 $\pm$ SD', alpha=.5, color='#4C4CA6')
axs[3,1].plot(data[:,0], averages_2[:,11] - stds_2[:,11], alpha=.5, color='#4C4CA6')

axs[3,1].axhline(num_exp[9], color='#FA8072', linewidth=3) # s_1-1_1 prot
axs[3,1].axhline(num_exp[9] + num_std[9], linestyle='--', color='#FA8072', linewidth=3)
axs[3,1].axhline(num_exp[9] - num_std[9], linestyle='--', color='#FA8072', linewidth=3)
axs[3,1].plot(data[:,0], averages_2[:,12], label='S_1-1_1', alpha=.8, color='#FA8072')
axs[3,1].plot(data[:,0], averages_2[:,12] + stds_2[:,12], label='S_1-1_1 $\pm$ SD', alpha=.5, color='#FA8072')
axs[3,1].plot(data[:,0], averages_2[:,12] - stds_2[:,12], alpha=.5, color='#FA8072')

axs[3,1].axhline(num_exp[10], color='#FFD700', linewidth=3) # s_1-1_2 prot
axs[3,1].axhline(num_exp[10] + num_std[10], linestyle='--', color='#FFD700', linewidth=3)
axs[3,1].axhline(num_exp[10] - num_std[10], linestyle='--', color='#FFD700', linewidth=3)
axs[3,1].plot(data[:,0], averages_2[:,13], label='S_1-1_2', alpha=.8, color='#FFD700')
axs[3,1].plot(data[:,0], averages_2[:,13] + stds_2[:,13], label='S_1-1_2 $\pm$ SD', alpha=.5, color='#FFD700')
axs[3,1].plot(data[:,0], averages_2[:,13] - stds_2[:,13], alpha=.5, color='#FFD700')

axs[3,1].axhline(num_exp[13], color='#7CFC00', linewidth=3) # s_1-2_1 prot
axs[3,1].axhline(num_exp[13] + num_std[13], linestyle='--', color='#7CFC00', linewidth=3)
axs[3,1].axhline(num_exp[13] - num_std[13], linestyle='--', color='#7CFC00', linewidth=3)
axs[3,1].plot(data[:,0], averages_2[:,14], label='S_1-2_1', alpha=.8, color='#7CFC00')
axs[3,1].plot(data[:,0], averages_2[:,14] + stds_2[:,14], label='S_1-2_1 $\pm$ SD', alpha=.5, color='#7CFC00')
axs[3,1].plot(data[:,0], averages_2[:,14] - stds_2[:,14], alpha=.5, color='#7CFC00')

axs[3,1].axhline(num_exp[14], color='#228B22', linewidth=3) # s_1-2_2 prot
axs[3,1].axhline(num_exp[14] + num_std[14], linestyle='--', color='#228B22', linewidth=3)
axs[3,1].axhline(num_exp[14] - num_std[14], linestyle='--', color='#228B22', linewidth=3)
axs[3,1].plot(data[:,0], averages_2[:,15], label='S_1-2_2', alpha=.8, color='#228B22')
axs[3,1].plot(data[:,0], averages_2[:,15] + stds_2[:,15], label='S_1-2_2 $\pm$ SD', alpha=.5, color='#228B22')
axs[3,1].plot(data[:,0], averages_2[:,15] - stds_2[:,15], alpha=.5, color='#228B22')

axs[3,1].set_facecolor('white')

for i in range(2):
    axs[3,i].set_xlabel(r'Time in hours', fontsize=12)
    for ax in axs[:,i]:
        # ax.set_xlim([0,num_exp_1[num_exp_1.shape[0]-1,0]])
        ax.set_xlim([0,1000])
        ax.set_ylim(0)
        ax.legend(loc=4, fontsize=6.15, handlelength=3).set_zorder(100)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# axs[0].set_title("Average counts over multiple trajectories", fontsize=12)

axs[0,0].set_ylabel('Active gene copies count', fontsize=12)
axs[1,0].set_ylabel('mRNA count', fontsize=12)
axs[2,0].set_ylabel('Protein count', fontsize=12)
axs[3,0].set_ylabel('Synaptic protein count', fontsize=12)

axs[0,0].set_ylim([0, 2.01])
axs[0,1].set_ylim([0, 2.01])

axs[1,0].set_ylim([0, 9.1])
axs[1,1].set_ylim([0, 9.1])

axs[2,0].set_ylim([0, 5.3e4])
axs[2,1].set_ylim([0, 5.3e4])

axs[3,0].set_ylim([0, 5.6e3])
axs[3,1].set_ylim([0, 5.6e3])




plt.tight_layout()

fig.patch.set_facecolor('none')
# plt.savefig('../../../../../../../conferences/cosyne2024/poster/img/stationary_num_MC_2_genes_new.png', dpi=300)
# plt.savefig('../../../../../../../../../oleg_windows/olegr/Desktop/COSYNE/stationary_num_MC_2_genes_new.png', dpi=300)

plt.show()

fig.clear()
