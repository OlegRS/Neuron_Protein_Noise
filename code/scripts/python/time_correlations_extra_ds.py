import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt


time_corr = np.genfromtxt("data/tc__extra_ds", delimiter=',')

# 0) soma__Gene: 1
# 1) soma__mRNA: 0.618683
# 2) d_1__mRNA: 1.99291
# 3) d_2__mRNA: 0.00484799
# 4) d_2-1__mRNA: 1.19178
# 5) d_2-2__mRNA: 1.19178

# 6) soma__Prot: 2445.18
# 7) d_1__Prot: 37450.6
# 8) s_1_1__Prot: 1984.72
# 9) s_1_2__Prot: 2153.13
# 10) d_2__Prot: 0.929753
# 11) d_2-1__Prot: 18475.2
# 12) s_d_2-1_1__Prot: 1420.86
# 13) s_d_2-1_2__Prot: 1382.05
# 14) d_2-2__Prot: 18437.7
# 15) s_d_2-2_1__Prot: 1516.03

# 16) s_d_2-2_2__Prot: 1510.53



d_tau = .01
x_max = 1500

# labels = ['s_1-2_2_prot__s_1-2_2_prot', 's_1-2_2_prot__s_1-2_1_prot',
#           's_1-2_2_prot__s_1-1_2_prot', 's_1-2_2_prot__s_1-1_1_prot',
#           's_1-2_2_prot__s_1_2_prot', 's_1-2_2_prot__s_1_1_prot',
#           's_1-2_2_prot__d_1-2_prot', 's_1-2_2_prot__d_1-1_prot','s_1-2_2_prot__d_1_prot',
#           's_1-2_2_prot__soma_prot'] 
plt.plot(time_corr[:int(1500/d_tau),0], time_corr[:int(1500/d_tau), 1:])#, label=labels)

# eigenvals = [12.4862,10.2006, 9.6934, 8.2924, 9.3113, 6.7968, 0.0040, 0.0227, 0.0244, 0.0620, 0.1925, 0.1112, 0.0615, 0.0432, 0.1667]
# for ev in eigenvals: 
#     plt.axvline(1/ev, color='grey')

plt.xlim([0,x_max])
plt.ylim([0,1.01])

plt.legend()
plt.xlabel(r'Delay ($\tau$) in hours')
plt.ylabel(r'Pearson Correlation Coefficient')

plt.show()
