import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
import math

# ############ PARAMETERS #############
x_lim_ = 700
# x_lim = x_lim_*10

# sim_run_count = 1000 # Number of files with Gillespie simulations
# mult_sim_run_count = sim_run_count
# time = 10
# #####################################

# num_file_names = ["../../data/soma_only/test"]
# data_file_names = ["../../data/soma_only/tests2/fast_gene_activation__fast_protein_decay_"]


fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

# # Finding index of the desired moment in time
# file_name = data_file_names[0] + '0'
# data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
# time_ind = 0

# for t in data[:, 0]:
#     if time>t:
#         time_ind+=1
#     else:
#         print('Time set to')
#         print('time=', data[time_ind, 1])
#         break
# print('time_ind=', time_ind)

# numbers_of_active_genes = []
# numbers_of_mRNAs = []
# numbers_of_proteins = []

# for file_id in range(sim_run_count):
#     file_name = data_file_names[0] + str(file_id)
#     print(file_name)
#     data = np.genfromtxt(file_name, delimiter=',')[1:x_lim+1, 1:]
#     numbers_of_active_genes.append(data[time_ind, 1])
#     numbers_of_mRNAs.append(data[time_ind, 2])
#     numbers_of_proteins.append(data[time_ind, 3])

# Loading numerical results
num = np.genfromtxt(num_file_names[0], delimiter=',')[:x_lim, :]
gene_mean = num[time_ind, 1]
gene_std = num[time_ind, 2]
mRNA_mean = num[time_ind, 3]
mRNA_std = num[time_ind, 4]
prot_mean = num[time_ind, 5]
prot_std = num[time_ind, 6]


x  = np.arange(0, 1000, .1)
x_max = 150
x_ = np.arange(0, x_max)
factorials = [math.factorial(x_[i]) for i in range(0, x_max)]
axs[0].hist(numbers_of_active_genes, density=True)#, bins = degrees, color='b', density=True, label="3-star")
axs[0].plot(x, 1/np.sqrt(2*np.pi*gene_std**2)*np.exp(-((x-gene_mean)/gene_std)**2/2))
# axs[0].plot(x_, gene_mean**x_/factorials*np.exp(-gene_mean))

axs[1].hist(numbers_of_mRNAs, density=True)#, bins = degrees, color='b', density=True, label="3-star")
axs[1].plot(x, 1/np.sqrt(2*np.pi*mRNA_std**2)*np.exp(-((x-mRNA_mean)/mRNA_std)**2/2))
# axs[1].plot(x_, mRNA_mean**x_/factorials*np.exp(-mRNA_mean))

axs[2].hist(numbers_of_proteins, density=True)#, bins = degrees, color='b', density=True, label="3-star")
axs[2].plot(x, 1/np.sqrt(2*np.pi*prot_std**2)*np.exp(-((x-prot_mean)/prot_std)**2/2))
# axs[2].plot(x_, prot_mean**x_/factorials*np.exp(-prot_mean))

for ax in axs:
    ax.grid()
    ax.set_xlim([0,x_lim_])
    ax.set_ylim(0)

axs[0].set_xlim([0,2])


axs[0].set_ylabel(r'Active gene counts', fontsize=20)
axs[1].set_ylabel(r'mRNA counts', fontsize=20)
axs[2].set_ylabel(r'Protein counts', fontsize=20)
axs[2].set_xlabel(r'Number', fontsize=20)

plt.show()
