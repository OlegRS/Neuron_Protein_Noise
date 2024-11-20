import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt


lambda_2 = transcription_rate = (3.*200/10000) * .001#*3600

mRNA_decay_rate = 1.2e-5#*3600
tau_2 = 1/mRNA_decay_rate

kappa = translation_rate = 0.021#*3600

protein_decay_rate = 1.21e-6#*3600
tau_3 = 1/protein_decay_rate

D_m = mRNA_diffusion_constant = 3.4e-3
D_p = protein_diffusion_constant = .24

mRNA_forward_trafficking_velocity = .5e-2 * 2
mRNA_backward_trafficking_velocity = .1e-2
v_m = mRNA_forward_trafficking_velocity - mRNA_backward_trafficking_velocity

protein_forward_trafficking_velocity = .5e-2/10
protein_backward_trafficking_velocity = 0
v_p = protein_forward_trafficking_velocity - protein_backward_trafficking_velocity


x_lim = 3000

lambda_m = (v_m - np.sqrt(v_m**2+4*D_m/tau_2))/(2*D_m)

n_exp = .5



fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

### mRNA concentrations
# xi_max = 10000
# n_comps = 200+1
# xi = np.linspace(0,xi_max, n_comps)

xi = np.linspace(0,x_lim, 1000)

R = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2)) * np.exp(lambda_m * xi)

# R_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp.csv", delimiter=',')
R_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_mRNA_double_mRNA_fwd_traf.csv", delimiter=',')

axs[0].plot(xi, R)
axs[0].plot(np.arange(0,len(R_discrete)/2, .5), R_discrete, label="R_discrete")


### Protein concentrations
Lambda_m = (v_p - np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)
Lambda_p = (v_p + np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)

P = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2))*kappa/((lambda_m-Lambda_m)*(Lambda_p-lambda_m)) * (np.exp(lambda_m * xi) + 1/Lambda_p*(lambda_m - v_p/D_p)*np.exp(Lambda_m * xi))

# P_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_prot.csv", delimiter=',')
# P_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_prot_fwd_traf.csv", delimiter=',')
P_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_prot_double_mRNA_fwd_traf.csv", delimiter=',')


axs[1].plot(xi, P)
axs[1].plot(np.arange(0,len(P_discrete)/2, .5), P_discrete/4.16, label="P_discrete") # /4.16

for ax in axs:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.legend()

axs[1].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=20)
axs[0].set_ylabel("mRNA concentration, " + r'$\mu m^{-3}$', fontsize=20)
axs[1].set_ylabel("Protein concentration, " + r'$\mu m^{-3}$', fontsize=20)

plt.tight_layout()
plt.show()
