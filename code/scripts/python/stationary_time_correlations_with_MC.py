import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os

# 0) t
# 1) soma__Gene: 1, -0

# 2) soma__mRNA: 1.7959, 1.34011
# 4) d_1__mRNA: 1.29813, 1.13935
# 8) d_1-1__mRNA: 0.952986, 0.97621
# 12) d_1-1__mRNA: 0.952986, 0.97621

# 3) soma__Prot: 7004, 1874.81
# 5) d_1__Prot: 2130.81, 668.535
# 9) d_1-1__Prot: 202.58, 146.263
# 13) d_1-1__Prot: 892.124, 568.153


# 0) soma__Gene: 1, 0.707107, 0.707107
# 1) soma__mRNA: 1.7959, 1.55908, 0.868134
# 2) d_1__mRNA: 1.29813, 1.22705, 0.945247
# 3) d_1-1__mRNA: 0.952986, 1.01336, 1.06335
# 4) d_1-2__mRNA: 0.952986, 1.01336, 1.06335
# 5) soma__Prot: 28966.8, 4699.66, 0.162243
# 6) d_1__Prot: 28522.7, 4693.1, 0.164539
# 7) s_1_1__Prot: 2087.86, 346.338, 0.165882
# 8) s_1_2__Prot: 2091.73, 346.976, 0.165879
# 9) d_1-1__Prot: 14643.7, 2837.9, 0.193797
# 10) s_d_1-1_1__Prot: 1396.25, 272.912, 0.195461
# 11) s_d_1-1_2__Prot: 990.123, 194.273, 0.196211
# 12) d_1-2__Prot: 14643.7, 2799.16, 0.191151
# 13) s_d_1-2_1__Prot: 1274.36, 245.971, 0.193016
# 14) s_d_1-2_2__Prot: 2487.71, 477.688, 0.192019

# 0) soma__Gene: 1, 2.2747e-05, 2.2747e-05
# 1) soma__mRNA: 1.7959, 1.34011, 0.746206
# 2) d_1__mRNA: 1.29813, 1.13935, 0.87769
# 3) d_1-1__mRNA: 0.952986, 0.97621, 1.02437
# 4) d_1-2__mRNA: 0.952986, 0.97621, 1.02437
# 5) soma__Prot: 28966.8, 3619.83, 0.124965
# 6) d_1__Prot: 28522.7, 3665.99, 0.128529
# 7) s_1_1__Prot: 2087.86, 271.931, 0.130244
# 8) s_1_2__Prot: 2091.73, 272.428, 0.13024
# 9) d_1-1__Prot: 14643.7, 2381.31, 0.162617
# 10) s_d_1-1_1__Prot: 1396.25, 229.817, 0.164596
# 11) s_d_1-1_2__Prot: 990.123, 163.852, 0.165486
# 12) d_1-2__Prot: 14643.7, 2343.53, 0.160037
# 13) s_d_1-2_1__Prot: 1274.36, 206.776, 0.162259
# 14) s_d_1-2_2__Prot: 2487.71, 400.703, 0.161073


main_dir = "../../data/gillespie/all_stationary_time_correlations/"


print("Loading numerical expectations and stds...")
time_corr = np.genfromtxt("data/tc_mc_full", delimiter=',')

################### Gillespie #####################
# 1 - Active genes

# 2 - Soma mRNA
# 4 - d_1 mRNA
# 6 - d_1_1 mRNAs
# 8 - d_1_2 mRNAs

# 3 - Soma proteins
# 5 - d_1 proteins
# 7 - d_1_1 proteins
# 9 - d_1_2 proteins

# 10 - s_1-1 proteins
# 11 - s_1-2 proteins
# 12 - s_1-1_1 proteins
# 13 - s_1-1_2 proteions
# 14 - s_1-2_1 proteins
# 15 - s_1-2_2 proteions

dir_count = len(os.listdir(main_dir))
file_count = 0
for i in range(dir_count):
    file_count += len(os.listdir(main_dir + "stationary_time_correlations_" + str(i)))
print("file_count=", file_count)

############ PARAMETERS #############
step = .01
max_sim_time = 7000
tau_max = 1500
n_points = int(tau_max/step)
time = max_sim_time - tau_max
tau_steps = int(tau_max/step)
# x_lim = n_points*step
sim_run_count = file_count
n_compartments = 15
#####################################

averages = np.zeros((n_points, n_compartments))
variances = np.zeros((n_points, n_compartments))
correlators = np.zeros((n_points, n_compartments))

for dir_name in os.listdir(main_dir): # Loop over directories
    dir_name = main_dir + dir_name
    for file_name in os.listdir(dir_name):
        file_name = dir_name + '/' + file_name
        print("Loading file: " + file_name)
        data = np.genfromtxt(file_name, delimiter=',')[int((time-tau_max)/step):int(time/step), 1:]
        for tau_ind in range(tau_steps):
            averages[tau_ind] += data[tau_steps-tau_ind-1, 1:]/sim_run_count
            variances[tau_ind] += data[tau_steps-tau_ind-1, 1:]**2/sim_run_count
            correlators[tau_ind] += data[tau_steps-tau_ind-1, data.shape[1]-1] * data[tau_steps-1, 1:]/sim_run_count


pccs = np.zeros((n_points, n_compartments))
for tau_ind in range(tau_steps):
    pccs[tau_ind] = (correlators[tau_ind] - averages[tau_ind, n_compartments-1]*averages[0])/np.sqrt((variances[tau_ind,n_compartments-1] - averages[tau_ind,n_compartments-1]**2)*(variances[0] - averages[0]**2))

plt.plot(pccs)
plt.plot(time_corr[:,0], time_corr[:,1:])
plt.show()
