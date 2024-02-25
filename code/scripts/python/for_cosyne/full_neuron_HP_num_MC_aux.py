fig, axs = plt.subplots(nrows=6, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

axs[0,0].set_ylabel(r'S_1-1', fontsize=15)
axs[1,0].set_ylabel(r'S_1-2', fontsize=15)
axs[2,0].set_ylabel(r'S_1-1_1', fontsize=15)
axs[3,0].set_ylabel(r'S_1-1_2', fontsize=15)
axs[4,0].set_ylabel(r'S_1-2_1', fontsize=15)
axs[5,0].set_ylabel(r'S_1-2_2', fontsize=15)

######################## 1st column ###########################
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6], label="Expectation", alpha=.7, color='#4CFFFF')
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6] + num_std_1[:,6], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6] - num_std_1[:,6], alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,0].plot(data[:,0]-5000+10, averages_1[:,10], label='Average active genes', alpha=.5, color='#4CFFFF')
axs[0,0].set_facecolor('lightgray')


axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7], label="Expectation", alpha=.7, color='#4C4CA6')
axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7] + num_std_1[:,7], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7] - num_std_1[:,7], alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,0].plot(data[:,0]-5000+10, averages_1[:,11], label='Average active genes', alpha=.5, color='#4C4CA6')

axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10], label="Expectation", alpha=.7, color='#FA8072')
axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10] + num_std_1[:,10], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FA8072')
axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10] - num_std_1[:,10], alpha=.7, linestyle='--', color='#FA8072')
axs[2,0].plot(data[:,0]-5000+10, averages_1[:,12], label='Average active genes', alpha=.5, color='#FA8072')

axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11], label="Expectation", alpha=.7, color='#FFD700')
axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11] + num_std_1[:,11], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FFD700')
axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11] - num_std_1[:,11], alpha=.7, linestyle='--', color='#FFD700')
axs[3,0].plot(data[:,0]-5000+10, averages_1[:,13], label='Average active genes', alpha=.5, color='#FFD700')

axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14], label="Expectation", alpha=.7, color='#7CFC00')
axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14] + num_std_1[:,14], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#7CFC00')
axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14] - num_std_1[:,14], alpha=.7, linestyle='--', color='#7CFC00')
axs[4,0].plot(data[:,0]-5000+10, averages_1[:,14], label='Average active genes', alpha=.5, color='#7CFC00')

axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15], label="Expectation", alpha=.7, color='#228B22')
axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15] + num_std_1[:,15], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#228B22')
axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15] - num_std_1[:,15], alpha=.7, linestyle='--', color='#228B22')
axs[5,0].plot(data[:,0]-5000+10, averages_1[:,15], label='Average active genes', alpha=.5, color='#228B22')

######################## 2nd column ###########################
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6], label="Expectation", alpha=.7, color='#4CFFFF')
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6] + num_std_2[:,6], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6] - num_std_2[:,6], alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,1].plot(data[:,0]-5000+10, averages_2[:,10], label='Average active genes', alpha=.5, color='#4CFFFF')

axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7], label="Expectation", alpha=.7, color='#4C4CA6')
axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7] + num_std_2[:,7], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7] - num_std_2[:,7], alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,1].plot(data[:,0]-5000+10, averages_2[:,11], label='Average active genes', alpha=.5, color='#4C4CA6')

axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10], label="Expectation", alpha=.7, color='#FA8072')
axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10] + num_std_2[:,10], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FA8072')
axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10] - num_std_2[:,10], alpha=.7, linestyle='--', color='#FA8072')
axs[2,1].plot(data[:,0]-5000+10, averages_2[:,12], label='Average active genes', alpha=.5, color='#FA8072')

axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11], label="Expectation", alpha=.7, color='#FFD700')
axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11] + num_std_2[:,11], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FFD700')
axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11] - num_std_2[:,11], alpha=.7, linestyle='--', color='#FFD700')
axs[3,1].plot(data[:,0]-5000+10, averages_2[:,13], label='Average active genes', alpha=.5, color='#FFD700')

axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14], label="Expectation", alpha=.7, color='#7CFC00')
axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14] + num_std_2[:,14], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#7CFC00')
axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14] - num_std_2[:,14], alpha=.7, linestyle='--', color='#7CFC00')
axs[4,1].plot(data[:,0]-5000+10, averages_2[:,14], label='Average active genes', alpha=.5, color='#7CFC00')

axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15], label="Expectation", alpha=.7, color='#228B22')
axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15] + num_std_2[:,15], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#228B22')
axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15] - num_std_2[:,15], alpha=.7, linestyle='--', color='#228B22')
axs[5,1].plot(data[:,0]-5000+10, averages_2[:,15], label='Average active genes', alpha=.5, color='#228B22')
axs[5,1].set_facecolor('lightgray')

######################## 3rd column ###########################
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6], label="Expectation", alpha=.7, color='#4CFFFF')
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6] + num_std_3[:,6], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6] - num_std_3[:,6], alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,2].plot(data[:,0]-5000+10, averages_3[:,10], label='Average active genes', alpha=.5, color='#4CFFFF')
axs[0,2].set_facecolor('lightgray')


axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7], label="Expectation", alpha=.7, color='#4C4CA6')
axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7] + num_std_3[:,7], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7] - num_std_3[:,7], alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,2].plot(data[:,0]-5000+10, averages_3[:,11], label='Average active genes', alpha=.5, color='#4C4CA6')

axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10], label="Expectation", alpha=.7, color='#FA8072')
axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10] + num_std_3[:,10], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FA8072')
axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10] - num_std_3[:,10], alpha=.7, linestyle='--', color='#FA8072')
axs[2,2].plot(data[:,0]-5000+10, averages_3[:,12], label='Average active genes', alpha=.5, color='#FA8072')

axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11], label="Expectation", alpha=.7, color='#FFD700')
axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11] + num_std_3[:,11], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FFD700')
axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11] - num_std_3[:,11], alpha=.7, linestyle='--', color='#FFD700')
axs[3,2].plot(data[:,0]-5000+10, averages_3[:,13], label='Average active genes', alpha=.5, color='#FFD700')

axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14], label="Expectation", alpha=.7, color='#7CFC00')
axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14] + num_std_3[:,14], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#7CFC00')
axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14] - num_std_3[:,14], alpha=.7, linestyle='--', color='#7CFC00')
axs[4,2].plot(data[:,0]-5000+10, averages_3[:,14], label='Average active genes', alpha=.5, color='#7CFC00')

axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15], label="Expectation", alpha=.7, color='#228B22')
axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15] + num_std_3[:,15], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#228B22')
axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15] - num_std_3[:,15], alpha=.7, linestyle='--', color='#228B22')
axs[5,2].plot(data[:,0]-5000+10, averages_3[:,15], label='Average active genes', alpha=.5, color='#228B22')


for i in range(3):
    axs[5,i].set_xlabel(r'Time in hours', fontsize=15)
    for ax in axs[:,i]:
        # ax.set_xlim([0,num_exp_1[num_exp_1.shape[0]-1,0]])
        ax.set_xlim([0,500])
        ax.set_ylim(0)
        # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        # ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=15)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
