fig, axs = plt.subplots(nrows=6, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

# ----- MC ------
# 10 - s_1-1 proteins
# 11 - s_1-2 proteins
# 12 - s_1-1_1 proteins
# 13 - s_1-1_2 proteions
# 14 - s_1-2_1 proteins
# 15 - s_1-2_2 proteions

# --------- Theory ---------
#   6) s_1_1__Prot: 198.77, 63.7868
#   7) s_1_2__Prot: 198.77, 63.7868
#   10) s_1_1-1__Prot: 18.8974, 14.256
#   11) s_1_1-2__Prot: 188.974, 136.465
#   14) s_1_2-1__Prot: 83.2206, 53.7035
#   15) s_1_2-2__Prot: 83.2206, 53.7035

axs[0,0].set_ylabel(r'S_1-1', fontsize=18)
axs[1,0].set_ylabel(r'S_1-2', fontsize=18)
axs[2,0].set_ylabel(r'S_1-1_1', fontsize=18)
axs[3,0].set_ylabel(r'S_1-1_2', fontsize=18)
axs[4,0].set_ylabel(r'S_1-2_1', fontsize=18)
axs[5,0].set_ylabel(r'S_1-2_2', fontsize=18)


axs[0,1].plot(num_exp[:,0], num_exp[:,6], label="Expectation", color='red', alpha=.7)
axs[0,1].plot(num_exp[:,0], num_exp[:,6] + num_std[:,6], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[0,1].plot(num_exp[:,0], num_exp[:,6] - num_std[:,6], color='red', alpha=.7, linestyle='--')
axs[0,1].plot(data[:,0]-5000+10, averages[:,10], label='Average active genes', color='red', alpha=.5)

axs[1,1].plot(num_exp[:,0], num_exp[:,7], label="Expectation", color='red', alpha=.7)
axs[1,1].plot(num_exp[:,0], num_exp[:,7] + num_std[:,7], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[1,1].plot(num_exp[:,0], num_exp[:,7] - num_std[:,7], color='red', alpha=.7, linestyle='--')
axs[1,1].plot(data[:,0]-5000+10, averages[:,11], label='Average active genes', color='red', alpha=.5)

axs[2,1].plot(num_exp[:,0], num_exp[:,10], label="Expectation", color='red', alpha=.7)
axs[2,1].plot(num_exp[:,0], num_exp[:,10] + num_std[:,10], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[2,1].plot(num_exp[:,0], num_exp[:,10] - num_std[:,10], color='red', alpha=.7, linestyle='--')
axs[2,1].plot(data[:,0]-5000+10, averages[:,12], label='Average active genes', color='red', alpha=.5)

axs[3,1].plot(num_exp[:,0], num_exp[:,11], label="Expectation", color='red', alpha=.7)
axs[3,1].plot(num_exp[:,0], num_exp[:,11] + num_std[:,11], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[3,1].plot(num_exp[:,0], num_exp[:,11] - num_std[:,11], color='red', alpha=.7, linestyle='--')
axs[3,1].plot(data[:,0]-5000+10, averages[:,13], label='Average active genes', color='red', alpha=.5)

axs[4,1].plot(num_exp[:,0], num_exp[:,14], label="Expectation", color='red', alpha=.7)
axs[4,1].plot(num_exp[:,0], num_exp[:,14] + num_std[:,14], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[4,1].plot(num_exp[:,0], num_exp[:,14] - num_std[:,14], color='red', alpha=.7, linestyle='--')
axs[4,1].plot(data[:,0]-5000+10, averages[:,14], label='Average active genes', color='red', alpha=.5)

axs[5,1].plot(num_exp[:,0], num_exp[:,15], label="Expectation", color='red', alpha=.7)
axs[5,1].plot(num_exp[:,0], num_exp[:,15] + num_std[:,15], label=r"Expectation $\pm$ std", color='red', alpha=.7, linestyle='--')
axs[5,1].plot(num_exp[:,0], num_exp[:,15] - num_std[:,15], color='red', alpha=.7, linestyle='--')
axs[5,1].plot(data[:,0]-5000+10, averages[:,15], label='Average active genes', color='red', alpha=.5)


# axs[1,0].set_ylabel(r'S_1-1_1', fontsize=18)
# axs[1,0].plot(num_exp[:,0], num_exp[:,7], label="Soma", color='red', alpha=.7)
# axs[1,0].plot(num_exp[:,0], num_exp[:,7] + num_std[:,7], label="Soma", color='red', alpha=.7, linestyle='--')
# axs[1,0].plot(num_exp[:,0], num_exp[:,7] - num_std[:,7], label="Soma", color='red', alpha=.7, linestyle='--')

# axs[2,0].set_ylabel(r'S_1-1_2', fontsize=18)
# axs[2,0].plot(num_exp[:,0], num_exp[:,11], label="Soma", color='blue', alpha=.7)
# axs[2,0].plot(num_exp[:,0], num_exp[:,11] + num_std[:,11], label="Soma", color='blue', alpha=.7, linestyle='--')
# axs[2,0].plot(num_exp[:,0], num_exp[:,11] - num_std[:,11], label="Soma", color='blue', alpha=.7, linestyle='--')


# axs[3,0].set_ylabel(r'S_1-2_1', fontsize=18)
# axs[3,0].plot(num_exp[:,0], num_exp[:,14], label="Soma", color='red', alpha=.7)
# axs[3,0].plot(num_exp[:,0], num_exp[:,14] + num_std[:,14], label="Soma", color='red', alpha=.7, linestyle='--')
# axs[3,0].plot(num_exp[:,0], num_exp[:,14] - num_std[:,14], label="Soma", color='red', alpha=.7, linestyle='--')

# axs[3,0].set_ylabel(r'S_1-2_2', fontsize=18)
# axs[3,0].plot(num_exp[:,0], num_exp[:,15], label="Soma", color='red', alpha=.7)
# axs[3,0].plot(num_exp[:,0], num_exp[:,15] + num_std[:,15], label="Soma", color='red', alpha=.7, linestyle='--')
# axs[3,0].plot(num_exp[:,0], num_exp[:,15] - num_std[:,15], label="Soma", color='red', alpha=.7, linestyle='--')




for i in range(3):
    axs[5,i].set_xlabel(r'Time in hours', fontsize=18)
    for ax in axs[:,i]:
        # ax.set_xlim([0,num_exp[num_exp.shape[0]-1,0]])
        ax.set_xlim([0,500])
        ax.set_ylim(0)
        # ax.legend(loc=4)

# axs[0].set_title("Average counts over multiple trajectories", fontsize=18)

plt.tight_layout()

# plt.savefig('../data/protein_numbers.png', dpi=300)


plt.show()

fig.clear()
