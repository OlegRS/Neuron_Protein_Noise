pccs = np.zeros((n_points, n_compartments))
for tau_ind in range(tau_steps):
    pccs[tau_ind] = (correlators[tau_ind] - averages[tau_ind, n_compartments-1]*averages[0])/np.sqrt((variances[tau_ind,n_compartments-1] - averages[tau_ind,n_compartments-1]**2)*(variances[0] - averages[0]**2))


taus = np.arange(0,step,tau_max)
plt.plot(time_corr[:150000,0], pccs)
plt.plot(time_corr[:,0], time_corr[:,1:])
plt.show()

np.save("./data/averages_5040", averages, allow_pickle=False)
np.save("./data/variances_5040", variances, allow_pickle=False)
np.save("./data/correlators_5040", correlators, allow_pickle=False)
np.save("./data/pccs_5040", pccs, allow_pickle=False)
