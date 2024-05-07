import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

# time_corr = np.genfromtxt("/home/oleg/sync/study/ulster/current/Neuron_Protein_Noise/code/bin/exe/tc", delimiter=',')

# time_corr = np.genfromtxt("/home/oleg/sync/study/ulster/current/Neuron_Protein_Noise/code/bin/exe/tc_full", delimiter=',')
# time_corr = np.genfromtxt("data/tc_full", delimiter=',')
# time_corr = np.genfromtxt("data/tc_full_different_synapses", delimiter=',')
# time_corr = np.genfromtxt("data/tc_full_syn_dec", delimiter=',')
# time_corr = np.genfromtxt("data/tc_single_always_active_gene", delimiter=',')
# time_corr = np.genfromtxt("data/tc_small_soma__syn_dec", delimiter=',')
# time_corr = np.genfromtxt("data/tc_small_soma__syn_dec__asymmetric_branching", delimiter=',')
time_corr = np.genfromtxt("data/tc_small_soma__syn_dec__asymmetric_branching_large_2nd_branch", delimiter=',')

# s_1_1__Prot: 1481.98
# s_1_2__Prot: 2254.64
# s_d_1-1_1__Prot: 785.879
# s_d_1-1_2__Prot: 1061.15
# s_d_1-2_1__Prot: 1221.7
# s_d_1-2_2__Prot: 951.099


# labels = ['s_1-2_prot__s_1-2_Prot', 's_1-2_prot__s_1-1_Prot', 's_1-2_prot__d_1_Prot'] 
# plt.plot(time_corr[:20000,0], time_corr[:20000, 1:], label=labels)



d_tau = .01
x_max = 1500

labels = ['s_1-2_2_prot__s_1-2_2_prot', 's_1-2_2_prot__s_1-2_1_prot',
          's_1-2_2_prot__s_1-1_2_prot', 's_1-2_2_prot__s_1-1_1_prot',
          's_1-2_2_prot__s_1_2_prot', 's_1-2_2_prot__s_1_1_prot',
          's_1-2_2_prot__d_1-2_prot', 's_1-2_2_prot__d_1-1_prot','s_1-2_2_prot__d_1_prot',
          's_1-2_2_prot__soma_prot'] 
plt.plot(time_corr[:int(1500/d_tau),0], time_corr[:int(1500/d_tau), 1:], label=labels)

# eigenvals = [12.4862,10.2006, 9.6934, 8.2924, 9.3113, 6.7968, 0.0040, 0.0227, 0.0244, 0.0620, 0.1925, 0.1112, 0.0615, 0.0432, 0.1667]
# for ev in eigenvals: 
#     plt.axvline(1/ev, color='grey')

plt.xlim([0,x_max])
plt.ylim([0,1.01])

plt.legend()
plt.xlabel(r'Delay ($\tau$) in hours')
plt.ylabel(r'Pearson Correlation Coefficient')

plt.show()
