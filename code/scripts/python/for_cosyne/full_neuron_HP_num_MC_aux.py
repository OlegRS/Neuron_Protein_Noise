fig, axs = plt.subplots(nrows=6, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.8))


axs[0,0].set_ylabel(r'$q_{11}$ (S_1-1)', fontsize=15)
axs[1,0].set_ylabel(r'$q_{12}$ (S_1-2)', fontsize=15)
axs[2,0].set_ylabel(r'$q_{21}$ (S_1-1_1)', fontsize=15)
axs[3,0].set_ylabel(r'$q_{22}$ (S_1-1_2)', fontsize=15)
axs[4,0].set_ylabel(r'$q_{31}$ (S_1-2_1)', fontsize=15)
axs[5,0].set_ylabel(r'$q_{32}$ (S_1-2_2)', fontsize=15)

######################## 1st column ###########################
num_exp_1[:,0] -= 10

axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6], label="Expectation", alpha=.7, color='#4CFFFF')
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6] + num_std_1[:,6], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,0].plot(num_exp_1[:,0], num_exp_1[:,6] - num_std_1[:,6], alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,0].plot(data[:,0]-5000, averages_1[:,10], label='Average active genes', alpha=.5, color='#4CFFFF')
axs[0,0].set_facecolor('lightgray')


axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7], label="Expectation", alpha=.7, color='#4C4CA6')
axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7] + num_std_1[:,7], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,0].plot(num_exp_1[:,0], num_exp_1[:,7] - num_std_1[:,7], alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,0].plot(data[:,0]-5000, averages_1[:,11], label='Average active genes', alpha=.5, color='#4C4CA6')
axs[1,0].set_facecolor('white')


axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10], label="Expectation", alpha=.7, color='#FA8072')
axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10] + num_std_1[:,10], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FA8072')
axs[2,0].plot(num_exp_1[:,0], num_exp_1[:,10] - num_std_1[:,10], alpha=.7, linestyle='--', color='#FA8072')
axs[2,0].plot(data[:,0]-5000, averages_1[:,12], label='Average active genes', alpha=.5, color='#FA8072')
axs[2,0].set_facecolor('white')

axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11], label="Expectation", alpha=.7, color='#FFD700')
axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11] + num_std_1[:,11], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FFD700')
axs[3,0].plot(num_exp_1[:,0], num_exp_1[:,11] - num_std_1[:,11], alpha=.7, linestyle='--', color='#FFD700')
axs[3,0].plot(data[:,0]-5000, averages_1[:,13], label='Average active genes', alpha=.5, color='#FFD700')
axs[3,0].set_facecolor('white')

axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14], label="Expectation", alpha=.7, color='#7CFC00')
axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14] + num_std_1[:,14], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#7CFC00')
axs[4,0].plot(num_exp_1[:,0], num_exp_1[:,14] - num_std_1[:,14], alpha=.7, linestyle='--', color='#7CFC00')
axs[4,0].plot(data[:,0]-5000, averages_1[:,14], label='Average active genes', alpha=.5, color='#7CFC00')
axs[4,0].set_facecolor('white')

axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15], label="Expectation", alpha=.7, color='#228B22')
axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15] + num_std_1[:,15], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#228B22')
axs[5,0].plot(num_exp_1[:,0], num_exp_1[:,15] - num_std_1[:,15], alpha=.7, linestyle='--', color='#228B22')
axs[5,0].plot(data[:,0]-5000, averages_1[:,15], label='Average active genes', alpha=.5, color='#228B22')
axs[5,0].set_facecolor('white')

######################## 2nd column ###########################
num_exp_2[:,0] -= 10

axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6], label="Expectation", alpha=.7, color='#4CFFFF')
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6] + num_std_2[:,6], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,1].plot(num_exp_2[:,0], num_exp_2[:,6] - num_std_2[:,6], alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,1].plot(data[:,0]-5000, averages_2[:,10], label='Average active genes', alpha=.5, color='#4CFFFF')
axs[0,1].set_facecolor('white')

axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7], label="Expectation", alpha=.7, color='#4C4CA6')
axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7] + num_std_2[:,7], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,1].plot(num_exp_2[:,0], num_exp_2[:,7] - num_std_2[:,7], alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,1].plot(data[:,0]-5000, averages_2[:,11], label='Average active genes', alpha=.5, color='#4C4CA6')
axs[1,1].set_facecolor('white')

axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10], label="Expectation", alpha=.7, color='#FA8072')
axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10] + num_std_2[:,10], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FA8072')
axs[2,1].plot(num_exp_2[:,0], num_exp_2[:,10] - num_std_2[:,10], alpha=.7, linestyle='--', color='#FA8072')
axs[2,1].plot(data[:,0]-5000, averages_2[:,12], label='Average active genes', alpha=.5, color='#FA8072')
axs[2,1].set_facecolor('white')

axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11], label="Expectation", alpha=.7, color='#FFD700')
axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11] + num_std_2[:,11], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FFD700')
axs[3,1].plot(num_exp_2[:,0], num_exp_2[:,11] - num_std_2[:,11], alpha=.7, linestyle='--', color='#FFD700')
axs[3,1].plot(data[:,0]-5000, averages_2[:,13], label='Average active genes', alpha=.5, color='#FFD700')
axs[3,1].set_facecolor('white')

axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14], label="Expectation", alpha=.7, color='#7CFC00')
axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14] + num_std_2[:,14], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#7CFC00')
axs[4,1].plot(num_exp_2[:,0], num_exp_2[:,14] - num_std_2[:,14], alpha=.7, linestyle='--', color='#7CFC00')
axs[4,1].plot(data[:,0]-5000, averages_2[:,14], label='Average active genes', alpha=.5, color='#7CFC00')
axs[4,1].set_facecolor('white')

axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15], label="Expectation", alpha=.7, color='#228B22')
axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15] + num_std_2[:,15], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#228B22')
axs[5,1].plot(num_exp_2[:,0], num_exp_2[:,15] - num_std_2[:,15], alpha=.7, linestyle='--', color='#228B22')
axs[5,1].plot(data[:,0]-5000, averages_2[:,15], label='Average active genes', alpha=.5, color='#228B22')
axs[5,1].set_facecolor('lightgray')

######################## 3rd column ###########################
num_exp_3[:,0] -= 10

axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6], label="Expectation", alpha=.7, color='#4CFFFF')
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6] + num_std_3[:,6], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,2].plot(num_exp_3[:,0], num_exp_3[:,6] - num_std_3[:,6], alpha=.7, linestyle='--', color='#4CFFFF')
axs[0,2].plot(data[:,0]-5000, averages_3[:,10], label='Average active genes', alpha=.5, color='#4CFFFF')
axs[0,2].set_facecolor('lightgray')


axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7], label="Expectation", alpha=.7, color='#4C4CA6')
axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7] + num_std_3[:,7], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,2].plot(num_exp_3[:,0], num_exp_3[:,7] - num_std_3[:,7], alpha=.7, linestyle='--', color='#4C4CA6')
axs[1,2].plot(data[:,0]-5000, averages_3[:,11], label='Average active genes', alpha=.5, color='#4C4CA6')
axs[1,2].set_facecolor('white')

axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10], label="Expectation", alpha=.7, color='#FA8072')
axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10] + num_std_3[:,10], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FA8072')
axs[2,2].plot(num_exp_3[:,0], num_exp_3[:,10] - num_std_3[:,10], alpha=.7, linestyle='--', color='#FA8072')
axs[2,2].plot(data[:,0]-5000, averages_3[:,12], label='Average active genes', alpha=.5, color='#FA8072')
axs[2,2].set_facecolor('white')

axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11], label="Expectation", alpha=.7, color='#FFD700')
axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11] + num_std_3[:,11], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#FFD700')
axs[3,2].plot(num_exp_3[:,0], num_exp_3[:,11] - num_std_3[:,11], alpha=.7, linestyle='--', color='#FFD700')
axs[3,2].plot(data[:,0]-5000, averages_3[:,13], label='Average active genes', alpha=.5, color='#FFD700')
axs[3,2].set_facecolor('white')

axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14], label="Expectation", alpha=.7, color='#7CFC00')
axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14] + num_std_3[:,14], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#7CFC00')
axs[4,2].plot(num_exp_3[:,0], num_exp_3[:,14] - num_std_3[:,14], alpha=.7, linestyle='--', color='#7CFC00')
axs[4,2].plot(data[:,0]-5000, averages_3[:,14], label='Average active genes', alpha=.5, color='#7CFC00')
axs[4,2].set_facecolor('white')

axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15], label="Expectation", alpha=.7, color='#228B22')
axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15] + num_std_3[:,15], label=r"Expectation $\pm$ std", alpha=.7, linestyle='--', color='#228B22')
axs[5,2].plot(num_exp_3[:,0], num_exp_3[:,15] - num_std_3[:,15], alpha=.7, linestyle='--', color='#228B22')
axs[5,2].plot(data[:,0]-5000, averages_3[:,15], label='Average active genes', alpha=.5, color='#228B22')
axs[5,2].set_facecolor('white')


for i in range(3):
    axs[5,i].set_xlabel(r'Time in hours', fontsize=15)
    for ax in axs[:,i]:
        # ax.set_xlim([0,num_exp_1[num_exp_1.shape[0]-1,0]])
        ax.set_xlim([-10,490])
        ax.set_ylim(0)
        # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        # ax.legend(loc=4)

# axs[0,0].set_title(r"$\eta_{1,1}\mapsto 10\eta_{1,1}$", fontsize=15)
# axs[0,1].set_title(r"$\eta_{3,2}\mapsto 10\eta_{3,2}$", fontsize=15)
# axs[0,2].set_title(r"$\eta_{1,1}\mapsto 10\eta_{1,1}$ and $\lambda_2\mapsto 2\lambda_2$", fontsize=15)

axs[0,0].set_title("10x increase in S_1_1 binding rate constant", fontsize=15)
axs[0,1].set_title("10x increase in S_1-2_2 binding rate constant", fontsize=15)
axs[0,2].set_title(r"10x in S_1_1 binding and 2x in transcription", fontsize=15)

axs[0,0].set_ylim([0,35000])
axs[0,1].set_ylim([0,6000])
axs[0,2].set_ylim([0,60000])
axs[1,0].set_ylim([0,6000])
axs[1,1].set_ylim([0,6000])
axs[1,2].set_ylim([0,6000])
axs[2,0].set_ylim([0,3200])
axs[2,1].set_ylim([0,3200])
axs[2,2].set_ylim([0,3250])
axs[3,0].set_ylim([0,3250])
axs[3,1].set_ylim([0,3250])
axs[3,2].set_ylim([0,3250])

axs[4,0].set_ylim([0,3250])
axs[4,1].set_ylim([0,3250])
axs[4,2].set_ylim([0,3250])



# ##### ZOOM ######
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# axins = zoomed_inset_axes(axs[0,0], 29, bbox_to_anchor=(.4, .47+.05, .6, .52+.05), bbox_transform=ax.transAxes, loc="lower left")
# # axins.plot(L, F)
# axins.axhline(0, color='black')
# # axins.grid()
# x1, x2, y1, y2 = 0, .01, -0.005, 0.005 # specify the limits
# axins.set_xlim(x1, x2) # apply the x-limits
# axins.set_ylim(y1, y2) # apply the y-limits




plt.tight_layout()

fig.patch.set_facecolor('none')
plt.savefig('../../../../../../../conferences/cosyne2024/poster/img/full_neuron_HP_num_MC.png', dpi=300)


plt.show()

fig.clear()

num_exp_1[:,0] += 10
num_exp_2[:,0] += 10
num_exp_3[:,0] += 10
