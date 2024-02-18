fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))


axs[0].plot(num[:,0], num[:,1], color='black', linewidth=2, zorder=1e4)
axs[0].plot(num[:,0], num[:,1] + num[:,2], linestyle='--', color='black', linewidth=2, zorder=1e4)
axs[0].plot(num[:,0], num[:,1] - num[:,2], linestyle='--', color='black', linewidth=2, zorder=1e4)


axs[0].set_ylabel(r'Active genes', fontsize=20)


axs[0].plot(averages[:,0], averages[:,1], label='Average active genes', color='grey', alpha=.5)
axs[0].plot(averages[:,0], averages[:,1] + variances[:,1], linestyle='--', color='grey', linewidth=2)
axs[0].plot(averages[:,0], averages[:,1] - variances[:,1], linestyle='--', color='grey', linewidth=2)


axs[1].plot(averages[:,0], averages[:,2], label='Soma', color='red', alpha=.5)
axs[1].plot(averages[:,0], averages[:,2] + variances[:,2], linestyle='--', color='red', linewidth=2)
axs[1].plot(averages[:,0], averages[:,2] - variances[:,2], linestyle='--', color='red', linewidth=2)
axs[1].plot(averages[:,0], averages[:,4], label="Dend_1", color='blue', alpha=.5)
axs[1].plot(averages[:,0], averages[:,4] + variances[:,4], linestyle='--', color='blue', linewidth=2)
axs[1].plot(averages[:,0], averages[:,4] - variances[:,4], linestyle='--', color='blue', linewidth=2)
axs[1].plot(averages[:,0], averages[:,6], label="Dend_1-1", color='orange', alpha=.5)
axs[1].plot(averages[:,0], averages[:,6] + variances[:,6], linestyle='--', color='orange', linewidth=2)
axs[1].plot(averages[:,0], averages[:,6] - variances[:,6], linestyle='--', color='orange', linewidth=2)
axs[1].plot(averages[:,0], averages[:,8], label="Dend_1-2", color='green', alpha=.5)
axs[1].plot(averages[:,0], averages[:,8] + variances[:,8], linestyle='--', color='green', linewidth=2)
axs[1].plot(averages[:,0], averages[:,8] - variances[:,8], linestyle='--', color='green', linewidth=2)
axs[1].set_ylabel(r'mRNAs', fontsize=20)



axs[1].plot(num[:,0], num[:,3], color='red', zorder=1e4, label='Soma_theor')
axs[1].plot(num[:,0], num[:,3] + num[:,4], linestyle='--', color='red')
axs[1].plot(num[:,0], num[:,3] - num[:,4], linestyle='--', color='red')
axs[1].plot(num[:,0], num[:,5], color='blue', zorder=1e4, label='Dend_1_theor')
axs[1].plot(num[:,0], num[:,5] + num[:,6], linestyle='--', color='blue')
axs[1].plot(num[:,0], num[:,5] - num[:,6], linestyle='--', color='blue')
axs[1].plot(num[:,0], num[:,7], color='orange', zorder=1e4, label='Dend_1-1_theor')
axs[1].plot(num[:,0], num[:,7] + num[:,8], linestyle='--', color='orange')
axs[1].plot(num[:,0], num[:,7] - num[:,8], linestyle='--', color='orange')
axs[1].plot(num[:,0], num[:,9], color='green', zorder=1e4, label='Dend_1-2_theor')
axs[1].plot(num[:,0], num[:,9] + num[:,10], linestyle='--', color='green')
axs[1].plot(num[:,0], num[:,9] - num[:,10], linestyle='--', color='green')




############# Bulk proteins #############
axs[2].plot(num[:,0], num[:,11], color='cyan', zorder=1e4, label='Soma_theor')
axs[2].plot(num[:,0], num[:,11] + num[:,12], linestyle='--', color='cyan')
axs[2].plot(num[:,0], num[:,11] - num[:,12], linestyle='--', color='cyan')
axs[2].plot(num[:,0], num[:,13], color='pink', zorder=1e4, label='Dend_1_theor')
axs[2].plot(num[:,0], num[:,13] + num[:,14], linestyle='--', color='pink')
axs[2].plot(num[:,0], num[:,13] - num[:,14], linestyle='--', color='pink')
axs[2].plot(num[:,0], num[:,19], color='brown', zorder=1e4, label='Dend_1_1_theor')
axs[2].plot(num[:,0], num[:,19] + num[:,20], linestyle='--', color='brown')
axs[2].plot(num[:,0], num[:,19] - num[:,20], linestyle='--', color='brown')

axs[2].plot(averages[:,0], averages[:,3], label="Soma", color='red', alpha=.7)
axs[2].plot(num[:,0], averages[:,3] + variances[:,3], linestyle='--', color='red')
axs[2].plot(num[:,0], averages[:,3] - variances[:,3], linestyle='--', color='red')

axs[2].plot(averages[:,0], averages[:,5], label="Dend_1", color='blue', alpha=.7)
axs[2].plot(averages[:,0], averages[:,5] + variances[:,5], linestyle='--', color='blue')
axs[2].plot(averages[:,0], averages[:,5] - variances[:,5], linestyle='--', color='blue')

axs[2].plot(averages[:,0], averages[:,7], label="Dend_1-1", color='orange', alpha=.7)
axs[2].plot(averages[:,0], averages[:,7] + variances[:,7], linestyle='--', color='orange')
axs[2].plot(averages[:,0], averages[:,7] - variances[:,7], linestyle='--', color='orange')

axs[2].plot(averages[:,0], averages[:,9], label="Dend_1-2", color='green', alpha=.7)
axs[2].plot(averages[:,0], averages[:,9] + variances[:,9], linestyle='--', color='green')
axs[2].plot(averages[:,0], averages[:,9] - variances[:,9], linestyle='--', color='green')

axs[2].set_ylabel(r'Bulk proteins', fontsize=20)



###################################################

axs[3].plot(num[:,0], num[:,15], color='maroon', zorder=1e4, label='Syn_1-1_theor')
axs[3].plot(num[:,0], num[:,15] + num[:,16], linestyle='--', color='maroon')
axs[3].plot(num[:,0], num[:,15] - num[:,16], linestyle='--', color='maroon')

axs[3].plot(num[:,0], num[:,21], color='red', zorder=1e4, label='Syn_1-1_1_theor')
axs[3].plot(num[:,0], num[:,21] + num[:,22], linestyle='--', color='red')
axs[3].plot(num[:,0], num[:,21] - num[:,22], linestyle='--', color='red')

axs[3].plot(num[:,0], num[:,27], color='red', zorder=1e4, label='Syn_1-2_1_theor')
axs[3].plot(num[:,0], num[:,27] + num[:,28], linestyle='--', color='red')
axs[3].plot(num[:,0], num[:,27] - num[:,28], linestyle='--', color='red')

axs[3].set_ylabel(r'Synaptic proteins', fontsize=18)



axs[3].plot(averages[:,0], averages[:,10], label="Synapse_1_1", color='maroon', alpha=.7)
axs[3].plot(num[:,0], averages[:,10] + variances[:,10], linestyle='--', color='maroon')
axs[3].plot(num[:,0], averages[:,10] - variances[:,10], linestyle='--', color='maroon')

# axs[3].plot(averages[:,0], averages[:,12], label="Synapse_1_2", alpha=.7)

axs[3].plot(averages[:,0], averages[:,12], label="Synapse_1-2_1", alpha=.5, color='red')
axs[3].plot(num[:,0], averages[:,12] + variances[:,12], linestyle='--', color='red')
axs[3].plot(num[:,0], averages[:,12] - variances[:,12], linestyle='--', color='red')

axs[3].plot(averages[:,0], averages[:,14], label="Synapse_1-2_1", alpha=.5, color='pink')
axs[3].plot(num[:,0], averages[:,14] + variances[:,14], linestyle='--', color='pink')
axs[3].plot(num[:,0], averages[:,14] - variances[:,14], linestyle='--', color='pink')

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
