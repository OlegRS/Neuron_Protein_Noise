import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

time_corr = np.genfromtxt("/home/oleg/sync/study/ulster/current/Neuron_Protein_Noise/code/bin/exe/tc", delimiter=',')

labels = ['s_1-2_prot__s_1-2_Prot', 's_1-2_prot__s_1-1_Prot', 's_1-2_prot__d_1_Prot'] 
plt.plot(time_corr[:20000,0], time_corr[:20000, 1:], label=labels)
plt.legend()
plt.xlabel(r'Delay ($\tau$) in hours')
plt.ylabel(r'Pearson Correlation Coefficient')

plt.show()
