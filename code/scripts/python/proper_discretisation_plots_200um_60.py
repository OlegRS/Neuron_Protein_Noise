import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import  matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

d1_mRNA_means = np.genfromtxt("data/simple_neuron_discretisation_mRNA_200um_60_d1.csv", delimiter=',')
d2_mRNA_means = np.genfromtxt("data/simple_neuron_discretisation_mRNA_200um_60_d2.csv", delimiter=',')
d3_mRNA_means = np.genfromtxt("data/simple_neuron_discretisation_mRNA_200um_60_d3.csv", delimiter=',')

d1_Prot_means = np.genfromtxt("data/simple_neuron_discretisation_Prot_200um_60_d1.csv", delimiter=',')
d2_Prot_means = np.genfromtxt("data/simple_neuron_discretisation_Prot_200um_60_d2.csv", delimiter=',')
d3_Prot_means = np.genfromtxt("data/simple_neuron_discretisation_Prot_200um_60_d3.csv", delimiter=',')

dendritic_length = 200


############### mRNA ##################
n_compartments = len(d1_mRNA_means)
axs[0,0].plot(np.linspace(0, dendritic_length, n_compartments), d1_mRNA_means/dendritic_length*n_compartments)
n_compartments = len(d2_mRNA_means)
axs[0,1].plot(np.linspace(0, dendritic_length, n_compartments), d2_mRNA_means/dendritic_length*n_compartments)
axs[0,2].plot(np.linspace(0, dendritic_length, n_compartments), d3_mRNA_means/dendritic_length*n_compartments)

axs[0,0].set_ylabel("mRNA_density", fontsize=15)

n_compartments = len(d1_Prot_means)
axs[1,0].plot(np.linspace(0, dendritic_length, n_compartments), d1_Prot_means/dendritic_length*n_compartments)
n_compartments = len(d2_Prot_means)
axs[1,1].plot(np.linspace(0, dendritic_length, n_compartments), d2_Prot_means/dendritic_length*n_compartments)
axs[1,2].plot(np.linspace(0, dendritic_length, n_compartments), d3_Prot_means/dendritic_length*n_compartments)

axs[1,0].set_ylabel("Protein_density", fontsize=15)

for i in range(3):
    axs[1,i].set_xlabel("Length", fontsize=15)

for i in range(3):
    axs[0,i].set_title("D_" + str(i+1), fontsize=15)



plt.tight_layout()
plt.show()
