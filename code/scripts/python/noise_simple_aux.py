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

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

############ PARAMETERS #############
x_lim = 1000
#####################################

num_file_names = ["../../data/soma_only/test"]
data_file_names = ["../../data/soma_only/tests/fast_gene_activation__fast_protein_decay"]
data_colour = ["black", "pink", "maroon", "lime", "cyan"]
num_colour = ["red", "red", "brown", "green", "blue"]
num_labels = ["Theor", "slow_genes_theor", "super_slow_genes", "slow_mRNAs_theor", "slow_prot_theor"]
data_labels = ["Sim", "slow_genes", "super_slow_genes", "slow_mRNAs", "slow_prot"]

fig_avrg, axs_avrg = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))
fig_var, axs_var = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))
fig_var_avrg, axs_var_avrg = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))


for i in range(np.array(data_file_names).shape[0]):
    num = np.genfromtxt(num_file_names[i], delimiter=',')[:x_lim*10, :]

    axs_avrg[0].plot(averages[:,0], averages[:,1], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_avrg[0].plot(num[:,0], num[:, 1], linestyle='-', color=num_colour[i], label=num_labels[i])
    axs_avrg[1].plot(averages[:,0], averages[:,2], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_avrg[1].plot(num[:,0], num[:,3], linestyle='-', color=num_colour[i], label=num_labels[i])
    axs_avrg[2].plot(averages[:,0], averages[:,3], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_avrg[2].plot(num[:,0], num[:,5], linestyle='-', color=num_colour[i], label=num_labels[i])
    
    axs_var[0].plot(variances[:,0], variances[:,1], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_var[0].plot(num[:,0], num[:, 2], linestyle='-', color=num_colour[i], label=num_labels[i])
    axs_var[1].plot(variances[:,0], variances[:,2], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_var[1].plot(num[:,0], num[:,4], linestyle='-', color=num_colour[i], label=num_labels[i])
    axs_var[2].plot(variances[:,0], variances[:,3], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_var[2].plot(num[:,0], num[:,6], linestyle='-', color=num_colour[i], label=num_labels[i])

    axs_var_avrg[0].plot(variances[:,0], variances[:,1]**2/averages[:,1], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_var_avrg[0].plot(num[:,0], num[:, 2]**2/num[:, 1], linestyle='-', color=num_colour[i], label=num_labels[i])
    axs_var_avrg[1].plot(variances[:,0], variances[:,2]**2/averages[:,2], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_var_avrg[1].plot(num[:,0], num[:,4]**2/num[:,3], linestyle='-', color=num_colour[i], label=num_labels[i])
    axs_var_avrg[2].plot(variances[:,0], variances[:,3]**2/averages[:,3], linestyle='--', color=data_colour[i], label=data_labels[i])
    axs_var_avrg[2].plot(num[:,0], num[:,6]**2/num[:,5], linestyle='-', color=num_colour[i], label=num_labels[i])



for i in range(0,3):
    axs_var[i].grid()
    axs_avrg[i].grid()
    axs_var_avrg[i].grid()

axs_avrg[0].set_ylabel(r'Number of active genes', fontsize=20)
axs_avrg[1].set_ylabel(r'mRNA counts', fontsize=20)
axs_avrg[2].set_ylabel(r'Protein counts', fontsize=20)
axs_avrg[2].set_xlabel(r'Time in hours', fontsize=20)
    
axs_var[0].set_ylabel(r'Number of active genes', fontsize=20)
axs_var[1].set_ylabel(r'mRNA counts', fontsize=20)
axs_var[2].set_ylabel(r'Protein counts', fontsize=20)
axs_var[2].set_xlabel(r'Time in hours', fontsize=20)
# axs_var[1].axhline(223.60679774997897)

axs_var_avrg[0].set_ylabel(r'Number of active genes', fontsize=20)
axs_var_avrg[1].set_ylabel(r'mRNA counts', fontsize=20)
axs_var_avrg[2].set_ylabel(r'Protein counts', fontsize=20)
axs_var_avrg[2].set_xlabel(r'Time in hours', fontsize=20)



for ax in axs_var:
    ax.set_xlim([0,x_lim])
    ax.set_ylim(0)
    ax.legend(loc=4)


axs_avrg[0].set_title(r"\bf{Averages}", fontsize=20)
axs_var[0].set_title(r"\bf{Standard deviations}", fontsize=20)
axs_var_avrg[0].set_title(r"\bf{Variance/Average ratios}", fontsize=20)

plt.tight_layout()
# plt.savefig("noise_simple.png", dpi=300)
plt.show()
