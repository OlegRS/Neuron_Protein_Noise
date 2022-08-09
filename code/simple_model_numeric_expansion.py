import numpy as np
from sympy import *

N = symbols('N') # Maximum number of active  gene replicas
l1p = symbols('lambda1^+') # Rate of gene activation 
l1m = symbols('lambda1^-') # Rate of gene deactivation
l2 = symbols('lambda2') # Transcription rate
tau1 = 1/(l1p+l1m)
tau2 = symbols('tau2')  # mRNA lifetime
tau3 = symbols('tau3')  # Protein lifetime

Pon = l1p/(l1p+l1m)

theta01 = symbols('theta01') # Rate of protein hopping from compartment 0 to compartment 1
theta10 = symbols('theta10') # Rate of protein hopping from compartment 1 to compartment 0
theta12 = symbols('theta12') # Rate of protein hopping from compartment 1 to compartment 2
theta21 = symbols('theta21') # Rate of protein hopping from compartment 2 to compartment 1
theta13 = symbols('theta13') # Rate of protein hopping from compartment 1 to compartment 3
theta31 = symbols('theta31') # Rate of protein hopping from compartment 3 to compartment 1

eta11 = symbols('eta11') # Rate of protein binding to 1st synapse of 1st dendrite
eta12 = symbols('eta12') # Rate of protein binding to 2nd synapse of 1st dendrite
eta21 = symbols('eta21') # Rate of protein binding to 1st synapse of 2nd dendrite
eta22 = symbols('eta22') # Rate of protein binding to 2nd synapse of 2nd dendrite
eta31 = symbols('eta31') # Rate of protein binding to 1st synapse of 3rd dendrite
eta32 = symbols('eta32') # Rate of protein binding to 2nd synapse of 3rd dendrite

gamma11 = symbols('gamma11') # Rate of protein unbinding from 1st synapse of 1st dendrite
gamma12 = symbols('gamma12') # Rate of protein unbinding from 2nd synapse of 1st dendrite
gamma21 = symbols('gamma21') # Rate of protein unbinding from 1st synapse of 2nd dendrite
gamma22 = symbols('gamma22') # Rate of protein unbinding from 2nd synapse of 2nd dendrite
gamma31 = symbols('gamma31') # Rate of protein unbinding from 1st synapse of 3rd dendrite
gamma32 = symbols('gamma32') # Rate of protein unbinding from 2nd synapse of 3rd dendrite

k0 = symbols('kappa0') # Translation rate in the soma
k1 = symbols('kappa1') # Translation rate in 1st dendritic compartment
k2 = symbols('kappa2') # Translation rate in 2nd dendritic compartment
k3 = symbols('kappa3') # Translation rate in 3rd dendritic compartment

nu01 = symbols('nu01') # mRNA hopping rate from the soma to the 1st dendrite
nu10 = symbols('nu10') # mRNA hopping rate from the 1st dendrite to the soma
nu12 = symbols('nu12') # mRNA hopping rate from the 1st dendritic compartment to the 2nd
nu21 = symbols('nu21') # mRNA hopping rate from the 2nd dendritic compartment to the 1st
nu13 = symbols('nu13') # mRNA hopping rate from the 1st dendritic compartment to the 3rd
nu31 = symbols('nu31') # mRNA hopping rate from the 3rd dendritic compartment to the 1st


######################################
# # Parameters
dend_length = 200. # um; length of dendritic segments

alpha_1 = (3*dend_length/10000)*0.001*3600 # /hour; mRNA transcription rate (0.001/s CaMKII Fonkeu)

beta_1 = 1.2e-5*3600 # /hour; mRNA degradation rate  (1.2e-5/s CaMKII Fonkeu)
alpha_2 = 0.021*3600 # /hour; protein synthesis rate  (0.021/s CaMKII Fonkeu)
beta_2 = 1.21e-6*3600 # /hour; protein degradation rate (1.21e-6/s CaMKII Fonkeu)
alpha_3 = 6e-1 # /hour; protein-synapse binding rate
beta_3 = 6. # + 2*randn(6) # /hour; protein-synapse unbinding rate
d_m = 3.4e-3*3600*(1/dend_length**2) # /hour; mRNA diffusion rate constant (3.4e-3 um2/s for CaMKII)
d_p =  0.24*3600*(1/dend_length**2) # /hour; protein diffusion rate constant (0.24 um2/s for CaMKII)
kf_m = 0.1*5e-2*3600*(1/dend_length) # /hour; mRNA forward trafficking rate  (4e-2 um/s for CaMKII net rate)
kb_m = 0.1*1e-2*3600*(1/dend_length) # /hour; mRNA backward trafficking rate (4e-2 um/s for CaMKII net rate)
kf_p = 0. # /hour; protein forward trafficking rate (0 for CaMKII)
kb_p = 0. # /hour; protein backward trafficking rate (0 for CaMKII)

# event_start = 14*24+100.
# event_duration = 30/60.
# psynth_event_rate = 1. # fold-increase in protein synthesis rate
# beta_3_event = beta_3[1]

######################################

N_ = 1
l1p_ = 1
l1m_ = 0 # This results to a single gene which is always active
l2_ = alpha_1
theta01_ = kf_p + d_p
theta10_ = kb_p + d_p
eta11_ = alpha_3
gamma11_ = beta_3

k0_ = alpha_2
k1_ = alpha_2
nu01_ = d_m + kf_m
nu10_ = d_m + kb_m

tau1_ = 1/(l1p_+l1m_)
tau2_ = 1/beta_1
tau3_ = 1/beta_2

Pon_ = l1p_/(l1p_+l1m_)


x = Matrix(symbols('x0:15'))

E = Matrix(symbols('E0:15'))   # Expectations (1st order)
GG  = symbols('GG')    # Gene-Gene (2nd order)
GM = symbols('GM0:4')  # Gene-mRNA correlations (2nd order)
GP = symbols('GP0:10') # Gene-Protein
MM = symbols('MM0:10') # mRNA-mRNA
MP = symbols('MP0:40') # mRNA-Protein
PP = symbols('PP0:55') # Protein-Protein

# Now we fill G_2 (2nd order) matrix

G_2 = Matrix(np.zeros((15,15)))

G_2[0,0] = GG

for i in range(1,5):
    G_2[0,i] = GM[i-1]
    G_2[i,0] = G_2[0,i]
    
for i in range(5,15):
    G_2[0,i] = GP[i-5]
    G_2[i,0] = G_2[0,i]

j_min, j_max = 1, 5
ind = 0
for i in np.arange(j_min, j_max):
    for j in np.arange(i, j_max):
        G_2[i,j] = MM[ind]
        G_2[j,i] = G_2[i,j]
        ind += 1

i_min, i_max = 1, 5
j_min, j_max = 5, 15
ind = 0
for i in np.arange(i_min, i_max):
    for j in np.arange(j_min, j_max):
        G_2[i,j] = MP[ind]
        G_2[j,i] = G_2[i,j]
        ind += 1

j_min, j_max = 5, 15
ind = 0
for i in np.arange(j_min, j_max):
    for j in np.arange(i, j_max):
        G_2[i,j] = PP[ind]
        G_2[j,i] = G_2[i,j]
        ind += 1

pretty_print(G_2)

G = 1 + (transpose(E)*x + 1/2*transpose(x)*G_2*x)[0,0]

############## DEFINING THE PDE ################
dGdx = Matrix(symbols('dGdx0:15'))
# G_ = symbols('G')

dGdt = l1p*(x[0]*N*G - (x[0]+1)*x[0]*dGdx[0]) - l1m*x[0]*dGdx[0] + l2*x[1]*(x[0]+1)*dGdx[0] - 1/tau2*(x[1]*dGdx[1]+x[2]*dGdx[2]+x[3]*dGdx[3]+x[4]*dGdx[4]) + nu01*(x[2]-x[1])*dGdx[1] + nu10*(x[1]-x[2])*dGdx[2] + nu12*(x[3]-x[2])*dGdx[2] + nu21*(x[2]-x[3])*dGdx[3] + nu13*(x[4]-x[2])*dGdx[2] + nu31*(x[2]-x[4])*dGdx[4] + k0*(x[1]+1)*x[5]*dGdx[1] + k1*(x[2]+1)*x[6]*dGdx[2] + k2*(x[3]+1)*x[7]*dGdx[3] + k3*(x[4]+1)*x[8]*dGdx[4] - 1/tau3*(x[5]*dGdx[5] + x[6]*dGdx[6] + x[7]*dGdx[7]) + theta01*(x[6]-x[5])*dGdx[5] + theta10*(x[5]-x[6])*dGdx[6] + theta12*(x[7]-x[6])*dGdx[6] + theta21*(x[6]-x[7])*dGdx[7] + theta13*(x[8]-x[6])*dGdx[6] + theta31*(x[6]-x[8])*dGdx[8] + eta11*(x[9]-x[6])*dGdx[6] + eta12*(x[10]-x[6])*dGdx[6] + eta21*(x[11]-x[7])*dGdx[7] + eta22*(x[12]-x[7])*dGdx[7] + eta31*(x[13]-x[8])*dGdx[8] + eta32*(x[14]-x[8])*dGdx[8] + gamma11*(x[6]-x[9])*dGdx[9] + gamma12*(x[6]-x[10])*dGdx[10] + gamma21*(x[7]-x[11])*dGdx[11] + gamma22*(x[7]-x[12])*dGdx[12] + gamma31*(x[8]-x[13])*dGdx[13] + gamma32*(x[8]-x[14])*dGdx[14] 


# A = list(symbols('A0:15')) # PDE coefficients
# A[0] = -collect(dGdt, dGdx[0], evaluate=False)[dGdx[0]]
# A[1] = -collect(dGdt, dGdx[1], evaluate=False)[dGdx[1]]
# A[2] = -collect(dGdt, dGdx[2], evaluate=False)[dGdx[2]]
# A[3] = -collect(dGdt, dGdx[3], evaluate=False)[dGdx[3]]
# A[4] = -collect(dGdt, dGdx[4], evaluate=False)[dGdx[4]]
# A[5] = -collect(dGdt, dGdx[5], evaluate=False)[dGdx[5]]
# A[6] = -collect(dGdt, dGdx[6], evaluate=False)[dGdx[6]]
# A[7] = -collect(dGdt, dGdx[7], evaluate=False)[dGdx[7]]
# A[8] = -collect(dGdt, dGdx[8], evaluate=False)[dGdx[8]]
# A[9] = -collect(dGdt, dGdx[9], evaluate=False)[dGdx[9]]
# A[10] = -collect(dGdt, dGdx[10], evaluate=False)[dGdx[10]]
# A[11] = -collect(dGdt, dGdx[11], evaluate=False)[dGdx[11]]
# A[12] = -collect(dGdt, dGdx[12], evaluate=False)[dGdx[12]]
# A[13] = -collect(dGdt, dGdx[13], evaluate=False)[dGdx[13]]
# A[14] = -collect(dGdt, dGdx[14], evaluate=False)[dGdx[14]]

# RHS = l1p*x[0]*N*G

dGdt = collect(dGdt, dGdx)

dGdt = dGdt.subs([(dGdx[i], diff(G, x)[i]) for i in range(0,15)])

# print('dGdt=')
# pretty_print(dGdt)

print("============ SOLVING FIRST ORDER... ===========")

X = Matrix(symbols('X0:15'))
for i in range(0,15):
    X[i] = (collect(expand(dGdt), x[i], evaluate=False)[x[i]]).subs([(x[j],0) for j in range(0, 15)])


# O1 = linsolve([X1, X2, X3, X4, X5, X6], G1, G2, G3, G4, G5, G6)

# G1_ = factor(O1.args[0][0])
# G2_ = factor(O1.args[0][1])
# G3_ = factor(O1.args[0][2])
# G4_ = factor(O1.args[0][3])
# G5_ = factor(O1.args[0][4])
# G6_ = factor(O1.args[0][5])

# print('G1=')
# pretty_print(G1_)
# print("------------------------------------------------")
# print('G2=')
# pretty_print(G2_)
# print("------------------------------------------------")
# print('G3=')
# pretty_print(G3_)
# print("------------------------------------------------")
# print('G4=')
# pretty_print(G4_)
# print("------------------------------------------------")
# print('G5=')
# pretty_print(G5_)
# print("------------------------------------------------")
# print('G6=')
# pretty_print(G6_)


