import sys
sys.path.append('../build')

import SGEN_Py as sg

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

#####################################################################################

Dendrite_length = 5000 #um
N_dendritic_segments = 100

soma = sg.Soma("soma", Dendrite_length/N_dendritic_segments)

dendritic_segments = [sg.Dendritic_segment(soma, "d_1-1", Dendrite_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    dendritic_segments.append(sg.Dendritic_segment(dendritic_segments[i-1], "d_1-" + str(i+1), Dendrite_length/N_dendritic_segments))
    
neuron = sg.Neuron(soma)

ae = sg.Analytic_engine(neuron)

print("Computing mRNA expectations...")
mRNA_expectations = np.array(ae.stationary_mRNA_expectations())
print("Computing protein expectations...")
prot_expectations = np.array(ae.stationary_protein_expectations())
# print("Computing gene-mRNA correlations...")
# gene_mRNA_covariances = np.array(ae.stationary_gene_mRNA_covariances())
# print("Computing mRNA-mRNA correlations...")
# mRNA_mRNA_covariances = np.array(ae.stationary_mRNA_mRNA_covariances())
# print("Computing gene-protein correlations...")
# gene_prot_covariances = np.array(ae.stationary_gene_protein_covariances())
# print("Computing mRNA-protein correlations...")
# mRNA_prot_covariances = np.array(ae.stationary_mRNA_protein_covariances())
# print("Computing protein-protein correlations...")
# prot_prot_covariances = np.array(ae.stationary_protein_protein_covariances())


soma = sg.Soma("soma", Dendrite_length/N_dendritic_segments)

dendritic_segments = [sg.Dendritic_segment(soma, "d_1-1", Dendrite_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    dendritic_segments.append(sg.Dendritic_segment(dendritic_segments[i-1], "d_1-" + str(i+1), Dendrite_length/N_dendritic_segments))
    
neuron = sg.Neuron(soma)


print("Computing expectations and correlations...")
# As for now, we need another analytic engine to compute again
# sg.Analytic_engine(neuron.refresh()).stationary_expectations_and_correlations()
sg.Analytic_engine(neuron).stationary_expectations_and_correlations()
# ae.stationary_expectations_and_correlations()

##############################################################################################
##############################################################################################

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

x_lim = Dendrite_length

lambda_m = (v_m - np.sqrt(v_m**2+4*D_m/tau_2))/(2*D_m)

n_exp = .5

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

### mRNA concentrations

xi = np.linspace(0,x_lim, 100000)

R = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2)) * np.exp(lambda_m * xi)
R_discrete = np.genfromtxt("mRNA_expectations.csv", delimiter=',')

axs[0,0].plot(xi, R, label="R_analytic")

axs[0,0].plot(np.arange(0,len(R_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), R_discrete * N_dendritic_segments/Dendrite_length, label="R_discrete")

axs[0,0].plot(np.arange(0,len(mRNA_expectations)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), mRNA_expectations * N_dendritic_segments/Dendrite_length, label="mRNA_concentrations")


### Protein concentrations
Lambda_m = (v_p - np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)
Lambda_p = (v_p + np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)

P = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2))*kappa/((lambda_m-Lambda_m)*(Lambda_p-lambda_m)*D_p) * (np.exp(lambda_m * xi) + 1/Lambda_p*(lambda_m - v_p/D_p)*np.exp(Lambda_m * xi))

P_discrete = 2*np.genfromtxt("protein_expectations.csv", delimiter=',')

axs[1,0].plot(xi, P, label="P_analytic")
axs[1,0].plot(np.arange(0,len(prot_expectations)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), P_discrete*N_dendritic_segments/Dendrite_length, label="P_discrete")

axs[1,0].plot(np.arange(0,len(prot_expectations)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), prot_expectations *N_dendritic_segments/Dendrite_length, label="Prot_concentration")


for ax in axs[:,0]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    # ax.set_yscale('log')
    ax.legend()

axs[1,0].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=20)
axs[0,0].set_ylabel("mRNA concentration, " + r'$\mu m^{-3}$', fontsize=20)
axs[1,0].set_ylabel("Protein concentration, " + r'$\mu m^{-3}$', fontsize=20)

### gene-mRNA correlations
lambda_m_tilde = (v_m - np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))/(2*D_m)

G2_nm = 2*lambda_2*n_exp*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))*np.exp(lambda_m_tilde*xi) + n_exp*R

G2_nm_discrete = np.genfromtxt("gene_mRNA_covariances.csv", delimiter=',')

axs[0,1].plot(xi, G2_nm, label="G2_nm")
axs[0,1].plot(np.arange(0,len(G2_nm_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_nm_discrete*N_dendritic_segments/Dendrite_length, label="G2_nm_discrete")

axs[0,1].set_ylabel(r"$G^{(2)}_{nm}$ density, " + r'$\mu m^{-3}$', fontsize=20)

### mRNA noise
G2_m2 = (2*lambda_2*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2))) + R[0])*R

axs[1,1].plot(xi, G2_m2, label=r"$G^{(2)}_{m^2}$")

mRNA_covariances = np.genfromtxt("mRNA_covariances.csv", delimiter=',')

# mRNA_STD = [np.sqrt(mRNA_covariances[i,i] - R_discrete[i]**2) for i in range(len(R_discrete))]

G2_m2_discrete = np.array([mRNA_covariances[i,i] - R_discrete[i] for i in range(len(R_discrete))])
axs[1,1].plot(np.arange(0,len(G2_m2_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_m2_discrete*(N_dendritic_segments/Dendrite_length)**2, label="G2_discrete_10")


axs[1,1].set_ylabel(r"$G^{(2)}_{m^2}$ density, " + r'$\mu m^{-3}$', fontsize=20)
axs[1,1].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=20)

for ax in axs[:,1]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()

### gene-prot correlations
G2_np_discrete = np.genfromtxt("gene_prot_covariances.csv", delimiter=',')

axs[0,2].plot(np.arange(0,len(G2_np_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_np_discrete*N_dendritic_segments/Dendrite_length, label="G2_np_discrete_10")

axs[0,2].set_ylabel(r"$G^{(2)}_{np}$ density, " + r'$\mu m^{-3}$', fontsize=20)

for ax in axs[:,2]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()

plt.tight_layout()
plt.show()

#################################

prot_expectations = np.genfromtxt("protein_expectations.csv", delimiter=',')
prot_correlations = np.genfromtxt("prot_prot_covariances.csv", delimiter=',')

prot_stds = np.zeros(prot_correlations.shape[0])
for i in range(prot_correlations.shape[0]):
    prot_stds[i] = np.sqrt(prot_correlations[i,i] - prot_expectations[i]*prot_expectations[i])

PCCs = np.zeros(prot_correlations.shape)
for i in range(PCCs.shape[0]):
    for j in range(i):
        PCCs[i,j] = PCCs[j,i] = (prot_correlations[i,j] - prot_expectations[i]*prot_expectations[j])/(prot_stds[i]*prot_stds[j])

plt.imshow(PCCs)

plt.show()
