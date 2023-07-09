import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

############ PARAMETERS #############
x_lim = 9999
sim_run_count = 1000 # Number of files with Gillespie simulations
mult_sim_run_count = sim_run_count
#####################################

fig2, axs2 = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

for i in range(np.array(data_file_names).shape[0]):

    axs2[0].plot(variances[:,0], variances[:,1], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs2[0].plot(num[:,0], num[:, 2], linestyle='-', color=num_colour[i], label=num_labels[i])

    axs2[1].plot(variances[:,0], variances[:,2], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs2[1].plot(num[:,0], num[:,4], linestyle='-', color=num_colour[i], label=num_labels[i])

    axs2[2].plot(variances[:,0], variances[:,3], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs2[2].plot(num[:,0], num[:,6], linestyle='-', color=num_colour[i], label=num_labels[i])


    
axs2[0].grid()
axs2[1].grid()
axs2[2].grid()
    

for ax in axs2:
    ax.set_xlim([0,x_lim])
    ax.set_ylim(0)
    ax.legend(loc=4)


plt.show()
