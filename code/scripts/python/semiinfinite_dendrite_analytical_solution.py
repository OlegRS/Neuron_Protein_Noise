import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

gene_activation_rate = 1/12/3600
gene_deactivation_rate = 1/12/3600

tau_1 = 1/(gene_activation_rate+gene_deactivation_rate)

lambda_2 = transcription_rate = (3.*200/10000) * .001#*3600

mRNA_decay_rate = 1.2e-5#*3600
tau_2 = 1/mRNA_decay_rate

kappa = translation_rate = 0.021#*3600

protein_decay_rate = 1.21e-6#*3600
tau_3 = 1/protein_decay_rate

D_m = mRNA_diffusion_constant = 3.4e-3
D_p = protein_diffusion_constant = .24

mRNA_forward_trafficking_velocity = .5e-2
mRNA_backward_trafficking_velocity = .1e-2
v_m = mRNA_forward_trafficking_velocity - mRNA_backward_trafficking_velocity

protein_forward_trafficking_velocity = 0 # .5e-2/10
protein_backward_trafficking_velocity = 0
v_p = protein_forward_trafficking_velocity - protein_backward_trafficking_velocity


x_lim = 5000

lambda_m = (v_m - np.sqrt(v_m**2+4*D_m/tau_2))/(2*D_m)

n_exp = gene_activation_rate*tau_1;

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

### mRNA concentrations
# xi_max = 10000
# n_comps = 200+1
# xi = np.linspace(0,xi_max, n_comps)

xi = np.linspace(0,x_lim, 401)#100000)

R = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2)) * np.exp(lambda_m * xi)

# R_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp.csv", delimiter=',')
# R_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_mRNA_double_mRNA_fwd_traf.csv", delimiter=',')
# R_discrete = np.genfromtxt("../../data/5000um_5000_comp_mRNA.csv", delimiter=',')
# R_discrete = np.genfromtxt("../../data/5000um_50_comp_mRNA.csv", delimiter=',')
# R_discrete = np.genfromtxt("../../data/5000um_100_comp_mRNA.csv", delimiter=',')
# R_discrete = np.genfromtxt("../../data/5000um_150_comp_mRNA.csv", delimiter=',')
R_discrete = np.genfromtxt("../../data/5000um_250_comp_mRNA.csv", delimiter=',')
R_discrete_ = np.genfromtxt("../../bin/exe/mRNA_expectations", delimiter=',')

axs[0,0].plot(xi, R, label="R_analytic")
axs[0,0].plot(np.arange(0,len(R_discrete)*5000/251-.001,5000/251), R_discrete*251/5000, label="R_discrete") # np.arange(0,len(R_discrete)/2, .5),
axs[0,0].plot(np.arange(0,len(R_discrete_)*5000/251-.001,5000/251), R_discrete_*251/5000, label="R_discrete_") # np.arange(0,len(R_discrete_)/2, .5),


### Protein concentrations
Lambda_m = (v_p - np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)
Lambda_p = (v_p + np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)

P = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2))*kappa/((lambda_m-Lambda_m)*(Lambda_p-lambda_m)*D_p) * (np.exp(lambda_m * xi) + 1/Lambda_p*(lambda_m - v_p/D_p)*np.exp(Lambda_m * xi))

# P_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_prot.csv", delimiter=',')
# P_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_prot_fwd_traf.csv", delimiter=',')
# P_discrete = 2*np.genfromtxt("../../data/5000um_10000_comp_prot_double_mRNA_fwd_traf.csv", delimiter=',')
# P_discrete = np.genfromtxt("../../data/5000um_5000_comp_prot.csv", delimiter=',')
# P_discrete = np.genfromtxt("../../data/5000um_50_comp_prot.csv", delimiter=',')
# P_discrete = np.genfromtxt("../../data/5000um_100_comp_prot.csv", delimiter=',')
P_discrete = np.genfromtxt("../../data/5000um_250_comp_prot.csv", delimiter=',')

axs[1,0].plot(xi, P, label="P_analytic")
axs[1,0].plot(np.arange(0,len(P_discrete)*5000/251-.001,5000/251), P_discrete*251/5000, label="P_discrete") # /4.16 np.arange(0,len(P_discrete)/2, .5),

for ax in axs[:,0]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()

axs[1,0].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=20)
axs[0,0].set_ylabel("mRNA concentration, " + r'$\mu m^{-3}$', fontsize=20)
axs[1,0].set_ylabel("Protein concentration, " + r'$\mu m^{-3}$', fontsize=20)


### gene-mRNA correlations
lambda_m_tilde = (v_m - np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))/(2*D_m)

G2_nm = 2*lambda_2*n_exp*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))*np.exp(lambda_m_tilde*xi) + n_exp*R

G2_nm_discrete = np.genfromtxt("../../data/5000um_250_comp_gene_mRNA_cov.csv", delimiter=',')
G2_nm_discrete_ = np.genfromtxt("../../bin/exe/gene_mRNA_covariances", delimiter=',')
G2_nm_discrete_150 = np.genfromtxt("../../data/5000um_150_comp_gene_mRNA_cov.csv", delimiter=',')
G2_nm_discrete_400 = np.genfromtxt("../../bin/exe/sem_gene_mRNA_covariances", delimiter=',')

G2_nm_discrete_250 = np.genfromtxt("../../bin/exe/sem_gene_mRNA_covariances_250", delimiter=',')

axs[0,1].plot(xi, G2_nm, label="G2_nm")
axs[0,1].plot(np.arange(0,len(G2_nm_discrete)*5000/251-.001,5000/251), G2_nm_discrete*251/5000, label="G2_nm_discrete_250")
# axs[0,1].plot(np.arange(0,len(G2_nm_discrete_)*5000/251-.001,5000/251), G2_nm_discrete_*251/5000, label="G2_nm_discrete_250_")
# axs[0,1].plot(np.arange(0,len(G2_nm_discrete_150)*5000/151-.001,5000/151), G2_nm_discrete_150*151/5000, label="G2_nm_discrete_150")
axs[0,1].plot(np.arange(0,len(G2_nm_discrete_400)*5000/401-.001,5000/401), G2_nm_discrete_400*401/5000, label="G2_nm_discrete_400")
# axs[0,1].plot(np.arange(0,len(G2_nm_discrete_250)*5000/251-.001,5000/251), G2_nm_discrete_250*251/5000, label="G2_nm_discrete_250")

axs[0,1].set_ylabel(r"$G^{(2)}_{nm}$ density, " + r'$\mu m^{-3}$', fontsize=20)

### mRNA noise
# G2_m2 = 2*lambda_2*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))*R + R**2
G2_m2 = (2*lambda_2*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2))) + R[0])*R
# G2_m2 = R**2

axs[1,1].plot(xi, G2_m2, label=r"$G^{(2)}_{m^2}$")

# axs[1,1].plot(xi, lambda_m*xi, label=r"$\lambda_m\xi$")
# axs[1,1].plot(xi, lambda_m_tilde*xi, label=r"$\tilde\lambda_m\xi$")

mRNA_STD = np.genfromtxt("../../data/5000um_250_comp_mRNA_mRNA_STD.csv", delimiter=',')
G2_m2_discrete = mRNA_STD**2 + R_discrete**2 - R_discrete
G2_m2_discrete_ = np.genfromtxt("../../bin/exe/sem_mRNA_covariances", delimiter=',')
# G2_m2_discrete_250 = np.genfromtxt("../../bin/exe/sem_mRNA_covariances_250", delimiter=',')
# G2_m2_discrete_400 = np.genfromtxt("../../bin/exe/sem_mRNA_covariances_400_always_active_gene", delimiter=',')
# axs[1,1].plot(np.arange(0,len(G2_m2_discrete)*5000/251-.001,5000/251), (G2_m2_discrete+R_discrete*251/5000)*(251/5000)**2, label="G2_discrete_250")
axs[1,1].plot(np.arange(0,len(G2_m2_discrete)*5000/251-.001,5000/251), G2_m2_discrete*(251/5000)**2, label="G2_discrete_250")

# axs[1,1].plot(np.arange(0,len(G2_m2_discrete_100)*5000/101-.001,5000/101), G2_m2_discrete_100*101/5000, label="G2_discrete_101")
R_discrete_ = np.genfromtxt("../../bin/exe/sem_mRNA_expectations", delimiter=',')
# axs[1,1].plot(np.arange(0,len(G2_m2_discrete_)*5000/401-.001,5000/401), (G2_m2_discrete_+R_discrete_*(401/5000))*(401/5000)**2, label="G2_discrete_400_")

axs[1,1].plot(np.arange(0,len(G2_m2_discrete_)*5000/401-.001,5000/401), G2_m2_discrete_*(401/5000)**2, label="G2_discrete_400_")

# axs[1,1].plot(np.arange(0,len(G2_m2_discrete_250)*5000/251-.001,5000/251), G2_m2_discrete_250*(251/5000)**2, label="G2_discrete_250")
# axs[1,1].plot(np.arange(0,len(G2_m2_discrete_400)*5000/401-.001,5000/401), G2_m2_discrete_400*(401/5000)**2, label="G2_discrete_400")


ratios = G2_m2_discrete_/G2_m2

axs[1,1].set_ylabel(r"$G^{(2)}_{m^2}$ density, " + r'$\mu m^{-3}$', fontsize=20)
axs[1,1].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=20)

for ax in axs[:,1]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()


plt.tight_layout()
plt.show()
