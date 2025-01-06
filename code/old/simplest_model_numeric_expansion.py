import numpy as np
from sympy import *

N = symbols('N') # Maximum number of active  gene replicas
l1p = symbols('lambda1^+') # Rate of gene activation 
l1m = symbols('lambda1^-') # Rate of gene deactivation
l2 = symbols('lambda2') # Transcription rate
theta01 = symbols('theta01') # Rate of protein hopping from compartment 0 to compartment 1
theta10 = symbols('theta10') # Rate of protein hopping from compartment 1 to compartment 0
eta11 = symbols('eta11') # Rate of protein binding to 1st synapse of 1st dendrite
gamma11 = symbols('gamma11') # Rate of protein unbinding from 1st synapse of 1st dendrite

k0 = symbols('kappa0') # Translation rate in the soma
k1 = symbols('kappa1') # Translation rate in 1st dendrite 
nu01 = symbols('nu01') # mRNA hopping rate from the soma to the 1st dendrite
nu10 = symbols('nu10') # mRNA hopping rate from the 1st dendrite to the soma

tau1 = 1/(l1p+l1m)
tau2 = symbols('tau2')  # mRNA lifetime
tau3 = symbols('tau3')  # Protein lifetime

Pon = l1p/(l1p+l1m)


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



x1 = symbols('x1')
x2 = symbols('x2')
x3 = symbols('x3')
x4 = symbols('x4')
x5 = symbols('x5')
x6 = symbols('x6')

G1 = symbols('G1')
G2 = symbols('G2')
G3 = symbols('G3')
G4 = symbols('G4')
G5 = symbols('G5')
G6 = symbols('G6')

G1 = symbols('G1')
G2 = symbols('G2')
G3 = symbols('G3')
G4 = symbols('G4')
G5 = symbols('G5')
G6 = symbols('G6')


G11 = symbols('G11')
G12 = symbols('G12')
G13 = symbols('G13')
G14 = symbols('G14')
G15 = symbols('G15')
G16 = symbols('G16')

G22 = symbols('G22')
G23 = symbols('G23')
G24 = symbols('G24')
G25 = symbols('G25')
G26 = symbols('G26')

G33 = symbols('G33')
G34 = symbols('G34')
G35 = symbols('G35')
G36 = symbols('G36')

G44 = symbols('G44')
G45 = symbols('G45')
G46 = symbols('G46')

G55 = symbols('G55')
G56 = symbols('G56')

G66 = symbols('G66')


G_1 = Matrix([[G1], [G2], [G3], [G4], [G5], [G6]])
G_2 = Matrix([[G11, G12, G13, G14, G15, G16],
              [G12, G22, G23, G24, G25, G26],
              [G13, G23, G33, G34, G35, G36],
              [G14, G24, G34, G44, G45, G46],
              [G15, G25, G35, G45, G55, G56],
              [G16, G26, G36, G46, G56, G66]
              ])

x = Matrix([[x1], [x2], [x3], [x4], [x5], [x6]])

G = 1 + (transpose(G_1)*x + 1/2*transpose(x)*G_2*x)[0,0]


a1 = l1p*x1**2 - l2*x2*x1 + x1/tau1 - l2*x2
a2 = -(k0*x2*x4 + k0*x4 + nu01*x3 - (1/tau2 + nu01)*x2)
a3 = -(k1*x3*x5 + k1*x5 + nu10*x2 - (nu10 + 1/tau2)*x3)
a4 = (1/tau3 + theta01)*x4 - theta01*x5
a5 = (1/tau3 + theta10 + eta11)*x5 - theta10*x4 - eta11*x6
a6 = gamma11*(x6-x5)

RHS = l1p*x1*N*G

LHS = a1*diff(G, x1) + a2*diff(G, x2) + a3*diff(G, x3) + a4*diff(G, x4) + a5*diff(G, x5) + a6*diff(G, x6) 




### pretty_print(expand(RHS - LHS))
### collect(expand(RHS - LHS), x1*x2, evaluate=False)[x1*x2]

print("-------------- FIRST ORDER ------------------\n")
print("_________________x1________________\n")
X1 = (collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
pretty_print(X1)
#print((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
print("\n")

print("_________________x2________________\n")
X2 = (collect(expand(RHS - LHS), x2, evaluate=False)[x2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
pretty_print(X2)
print("\n")

print("_________________x3________________\n")
X3 = (collect(expand(RHS - LHS), x3, evaluate=False)[x3]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
pretty_print(X3)
print("\n")

print("_________________x4________________\n")
X4 = (collect(expand(RHS - LHS), x4, evaluate=False)[x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
pretty_print(X4)
print("\n")

print("_________________x5________________\n")
X5 = (collect(expand(RHS - LHS), x5, evaluate=False)[x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
pretty_print(X5)
print("\n")

print("_________________x6________________\n")
X6 = (collect(expand(RHS - LHS), x6, evaluate=False)[x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
pretty_print(X6)
print("\n")


print("============ SOLVING FIRST ORDER... ===========")

O1 = linsolve([X1, X2, X3, X4, X5, X6], G1, G2, G3, G4, G5, G6)

G1_ = factor(O1.args[0][0])
G2_ = factor(O1.args[0][1])
G3_ = factor(O1.args[0][2])
G4_ = factor(O1.args[0][3])
G5_ = factor(O1.args[0][4])
G6_ = factor(O1.args[0][5])

print('G1=')
pretty_print(G1_)
print("------------------------------------------------")
print('G2=')
pretty_print(G2_)
print("------------------------------------------------")
print('G3=')
pretty_print(G3_)
print("------------------------------------------------")
print('G4=')
pretty_print(G4_)
print("------------------------------------------------")
print('G5=')
pretty_print(G5_)
print("------------------------------------------------")
print('G6=')
pretty_print(G6_)


print("\n========= FIRST ORDER NUMERICAL SOLUTIONS: =========")
G1_num = G1_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
print("G1=<n>=", G1_num)
G2_num = G2_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
print("G2=<m0>=", G2_num)
G3_num = G3_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
print("G3=<m1>=", G3_num)
G4_num = G4_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
print("G4=<l0>=", G4_num)
G5_num = G5_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
print("G5=<l1>=", G5_num)
G6_num = G6_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
print("G6=<q11>=", G6_num)


##################### SECOND ORDER #########################

print("-------------- SECOND ORDER ------------------\n")

print("_________________x1^2________________\n")
X11 = (collect(expand(RHS - LHS), x1**2, evaluate=False)[x1**2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X11 = collect(X11, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X11)
print("\n")

print("_________________x1*x2________________\n")
X12 = (collect(expand(RHS - LHS), x1*x2, evaluate=False)[x1*x2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X12 = collect(X12, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X12)
print("\n")

print("_________________x1*x3________________\n")
X13 = (collect(expand(RHS - LHS), x1*x3, evaluate=False)[x1*x3]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X13 = collect(X13, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X13)
print("\n")

print("_________________x1*x4________________\n")
X14 = (collect(expand(RHS - LHS), x1*x4, evaluate=False)[x1*x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X14 = collect(X14, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X14)
print("\n")

print("_________________x1*x5________________\n")
X15 = (collect(expand(RHS - LHS), x1*x5, evaluate=False)[x1*x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X15 = collect(X15, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X15)
print("\n")

print("_________________x1*x6________________\n")
X16 = (collect(expand(RHS - LHS), x1*x6, evaluate=False)[x1*x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X16 = collect(X16, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X16)
print("\n")


print("_________________x2^2________________\n")
X22 = (collect(expand(RHS - LHS), x2**2, evaluate=False)[x2**2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X22 = collect(X22, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X22)
print("\n")

print("_________________x2*x3________________\n")
X23 = (collect(expand(RHS - LHS), x2*x3, evaluate=False)[x2*x3]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X23 = collect(X23, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X23)
print("\n")

print("_________________x2*x4________________\n")
X24 = (collect(expand(RHS - LHS), x2*x4, evaluate=False)[x2*x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X24 = collect(X24, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X24)
print("\n")

print("_________________x2*x5________________\n")
X25 = (collect(expand(RHS - LHS), x2*x5, evaluate=False)[x2*x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X25 = collect(X25, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X25)
print("\n")

print("_________________x2*x6________________\n")
X26 = (collect(expand(RHS - LHS), x2*x6, evaluate=False)[x2*x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X26 = collect(X26, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X26)
print("\n")


print("_________________x3^2________________\n")
X33 = (collect(expand(RHS - LHS), x3**2, evaluate=False)[x3**2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X33 = collect(X33, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X33)
print("\n")

print("_________________x3*x4________________\n")
X34 = (collect(expand(RHS - LHS), x3*x4, evaluate=False)[x3*x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X34 = collect(X34, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X34)
print("\n")

print("_________________x3*x5________________\n")
X35 = (collect(expand(RHS - LHS), x3*x5, evaluate=False)[x3*x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X35 = collect(X35, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X35)
print("\n")

print("_________________x3*x6________________\n")
X36 = (collect(expand(RHS - LHS), x3*x6, evaluate=False)[x3*x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X36 = collect(X36, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X36)
print("\n")


print("_________________x4^2________________\n")
X44 = (collect(expand(RHS - LHS), x4**2, evaluate=False)[x4**2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X44 = collect(X44, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X44)
print("\n")

print("_________________x4*x5________________\n")
X45 = (collect(expand(RHS - LHS), x4*x5, evaluate=False)[x4*x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X45 = collect(X45, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X45)
print("\n")

print("_________________x4*x6________________\n")
X46 = (collect(expand(RHS - LHS), x4*x6, evaluate=False)[x4*x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X46 = collect(X46, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X46)
print("\n")


print("_________________x5^2________________\n")
X55 = (collect(expand(RHS - LHS), x5**2, evaluate=False)[x5**2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X55 = collect(X55, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X55)
print("\n")

print("_________________x5*x6________________\n")
X56 = (collect(expand(RHS - LHS), x5*x6, evaluate=False)[x5*x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X56 = collect(X56, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X56)
print("\n")


print("_________________x6^2________________\n")
X66 = (collect(expand(RHS - LHS), x6**2, evaluate=False)[x6**2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)])
X66 = collect(X66, [G1, G2, G3, G4, G5, G6, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66]) 
pretty_print(X66)
print("\n")


# G_2_11 = collect(X11, G11, evaluate=False)[G11]
# G_2_12 = 0 #collect(X11, G12, evaluate=False)[G12]
# G_2_13 = 0 #collect(X11, G13, evaluate=False)[G13]
# G_2_14 = 0 #collect(X11, G14, evaluate=False)[G14]
# G_2_15 = 0 #collect(X11, G15, evaluate=False)[G15]
# G_2_16 = 0 #collect(X11, G16, evaluate=False)[G16]

# G_2_21 = G_2_12
# G_2_22 = collect(X11, G22, evaluate=False)[G22]
# G_2_23 = collect(X11, G23, evaluate=False)[G23]
# G_2_24 = collect(X11, G24, evaluate=False)[G24]
# G_2_25 = collect(X11, G25, evaluate=False)[G25]
# G_2_26 = collect(X11, G26, evaluate=False)[G26]

# G_2_31 = G_2_13
# G_2_32 = G_2_23
# G_2_33 = collect(X11, G33, evaluate=False)[G33]
# G_2_34 = collect(X11, G34, evaluate=False)[G34]
# G_2_35 = collect(X11, G35, evaluate=False)[G35]
# G_2_36 = collect(X11, G36, evaluate=False)[G36]

# G_2_41 = G_2_14
# G_2_42 = G_2_24
# G_2_43 = G_2_34
# G_2_44 = collect(X11, G44, evaluate=False)[G44]
# G_2_45 = collect(X11, G45, evaluate=False)[G45]
# G_2_46 = collect(X11, G46, evaluate=False)[G46]

# G_2_51 = G_2_15
# G_2_52 = G_2_25
# G_2_53 = G_2_35
# G_2_54 = G_2_45
# G_2_55 = collect(X11, G55, evaluate=False)[G55]
# G_2_56 = collect(X11, G56, evaluate=False)[G56]

# G_2_61 = G_2_16
# G_2_62 = G_2_26
# G_2_63 = G_2_36
# G_2_64 = G_2_46
# G_2_65 = G_2_56
# G_2_66 = collect(X11, G66, evaluate=False)[G66]



# G_2 = Matrix([[G_2_11, G_2_12, G_2_13, G_2_14, G_2_15, G_2_16],
#               [G_2_12, G_2_22, G_2_23, G_2_24, G_2_25, G_2_26],
#               [G_2_13, G_2_23, G_2_33, G_2_34, G_2_35, G_2_36],
#               [G_2_14, G_2_24, G_2_34, G_2_44, G_2_45, G_2_46],
#               [G_2_15, G_2_25, G_2_35, G_2_45, G_2_55, G_2_56],
#               [G_2_16, G_2_26, G_2_36, G_2_46, G_2_56, G_2_66]
#               ])


print("\n============ SOLVING EQUATION 1 FOR ACTIVE GENE VARIANCE... ==============")
O2_1_2_3 = linsolve([X11], G11)

G11_ = O2_1_2_3.args[0][0]

G11__ = G11_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G11_num = G11__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num)])
print("G11=", G11__, '=', G11_num)
print("Active gene replicas variance:")
print("var(n)=<n^2> - <n>^2 = G11 - G1*(G1-1)", G11_num - G1_num*(G1_num-1))

print("\n============ SOLVING EQUATIONS 2,3 FOR GENE-mRNA CORRELATIONS... ==============")
O2_1_2_3 = linsolve([X12, X13], G12, G13)
G12_ = O2_1_2_3.args[0][0]
G13_ = O2_1_2_3.args[0][1]

G12__ = G12_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G12_num = G12__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num)])

G13__ = G13_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G13_num = G13__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num)])
print("G12=", G12__, '=', G12_num)
print("G13=", G13__, '=', G13_num)

print("\n============ SOLVING EQUATIONS 4,5,6 FOR GENE-PROTEIN CORRELATIONS... ==============")
O2_4_5_6 = linsolve([X14, X15, X16], G14, G15, G16)

G14_ = O2_4_5_6.args[0][0]
G15_ = O2_4_5_6.args[0][1]
G16_ = O2_4_5_6.args[0][2]

G14__ = G14_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G14_num = G14__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num)])

G15__ = G15_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G15_num = G15__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num)])

G16__ = G16_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G16_num = G16__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num)])

print("G14=", G14__, '=', G14_num)
print("G15=", G15__, '=', G15_num)
print("G16=", G16__, '=', G16_num)


print("\n======== SOLVING EQUATIONS 7,8,12 FOR mRNA-mRNA CORRELATIONS AND VARIANCES... =========")
O2_7_8_12 = linsolve([X22, X23, X33], G22, G23, G33)

G22_ = O2_7_8_12.args[0][0]
G23_ = O2_7_8_12.args[0][1]
G33_ = O2_7_8_12.args[0][2]

G22__ = G22_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G22_num = G22__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num)])


G23__ = G23_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G23_num = G23__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num)])


G33__ = G33_.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_)])
G33_num = G33__.subs([(G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num)])

print("G22=", G22__, '=', G22_num)
print("G23=", G23__, '=', G23_num)
print("G33=", G33__, '=', G33_num)

print("mRNA variances:")
print("Soma: var(m_0)=<m_0^2> - <m_0>^2 = G22 - G2*(G2-1)", G22_num - G2_num*(G2_num-1))
print("Dendrite: var(m_1)=<m_1^2> - <m_1>^2 = G33 - G3*(G3-1)", G33_num - G3_num*(G3_num-1))


print("\n============ SOLVING EQUATIONS 9,10,11,13,14,15 FOR mRNA-PROTEIN CORRELATIONS... ==============")
# O2_9_10_11_13_14_15 = linsolve([X24, X25, X26, X34, X35, X36], G24, G25, G26, G34, G35, G36)

# G24 = O2_9_10_11_13_14_15.args[0][0]
# G25 = O2_9_10_11_13_14_15.args[0][1]
# G26 = O2_9_10_11_13_14_15.args[0][2]
# G34 = O2_9_10_11_13_14_15.args[0][3]
# G35 = O2_9_10_11_13_14_15.args[0][4]
# G36 = O2_9_10_11_13_14_15.args[0][5]

############# INVERTING MATRIX ############ 
M11 = collect(X24, G24, evaluate=False)[G24] if G24 in collect(X24, G24, evaluate=False) else 0
M12 = collect(X24, G25, evaluate=False)[G25] if G25 in collect(X24, G25, evaluate=False) else 0
M13 = collect(X24, G26, evaluate=False)[G26] if G26 in collect(X24, G26, evaluate=False) else 0
M14 = collect(X24, G34, evaluate=False)[G34] if G34 in collect(X24, G34, evaluate=False) else 0
M15 = collect(X24, G35, evaluate=False)[G35] if G35 in collect(X24, G35, evaluate=False) else 0
M16 = collect(X24, G36, evaluate=False)[G36] if G36 in collect(X24, G36, evaluate=False) else 0

M21 = collect(X25, G24, evaluate=False)[G24] if G24 in collect(X25, G24, evaluate=False) else 0
M22 = collect(X25, G25, evaluate=False)[G25] if G25 in collect(X25, G25, evaluate=False) else 0
M23 = collect(X25, G26, evaluate=False)[G26] if G26 in collect(X25, G26, evaluate=False) else 0
M24 = collect(X25, G34, evaluate=False)[G34] if G34 in collect(X25, G34, evaluate=False) else 0
M25 = collect(X25, G35, evaluate=False)[G35] if G35 in collect(X25, G35, evaluate=False) else 0
M26 = collect(X25, G36, evaluate=False)[G36] if G36 in collect(X25, G36, evaluate=False) else 0

M31 = collect(X26, G24, evaluate=False)[G24] if G24 in collect(X26, G24, evaluate=False) else 0
M32 = collect(X26, G25, evaluate=False)[G25] if G25 in collect(X26, G25, evaluate=False) else 0
M33 = collect(X26, G26, evaluate=False)[G26] if G26 in collect(X26, G26, evaluate=False) else 0
M34 = collect(X26, G34, evaluate=False)[G34] if G34 in collect(X26, G34, evaluate=False) else 0
M35 = collect(X26, G35, evaluate=False)[G35] if G35 in collect(X26, G35, evaluate=False) else 0
M36 = collect(X26, G36, evaluate=False)[G36] if G36 in collect(X26, G36, evaluate=False) else 0

M41 = collect(X34, G24, evaluate=False)[G24] if G24 in collect(X34, G24, evaluate=False) else 0
M42 = collect(X34, G25, evaluate=False)[G25] if G25 in collect(X34, G25, evaluate=False) else 0
M43 = collect(X34, G26, evaluate=False)[G26] if G26 in collect(X34, G26, evaluate=False) else 0
M44 = collect(X34, G34, evaluate=False)[G34] if G34 in collect(X34, G34, evaluate=False) else 0
M45 = collect(X34, G35, evaluate=False)[G35] if G35 in collect(X34, G35, evaluate=False) else 0
M46 = collect(X34, G36, evaluate=False)[G36] if G36 in collect(X34, G36, evaluate=False) else 0

M51 = collect(X35, G24, evaluate=False)[G24] if G24 in collect(X35, G24, evaluate=False) else 0
M52 = collect(X35, G25, evaluate=False)[G25] if G25 in collect(X35, G25, evaluate=False) else 0
M53 = collect(X35, G26, evaluate=False)[G26] if G26 in collect(X35, G26, evaluate=False) else 0
M54 = collect(X35, G34, evaluate=False)[G34] if G34 in collect(X35, G34, evaluate=False) else 0
M55 = collect(X35, G35, evaluate=False)[G35] if G35 in collect(X35, G35, evaluate=False) else 0
M56 = collect(X35, G36, evaluate=False)[G36] if G36 in collect(X35, G36, evaluate=False) else 0

M61 = collect(X36, G24, evaluate=False)[G24] if G24 in collect(X36, G24, evaluate=False) else 0
M62 = collect(X36, G25, evaluate=False)[G25] if G25 in collect(X36, G25, evaluate=False) else 0
M63 = collect(X36, G26, evaluate=False)[G26] if G26 in collect(X36, G26, evaluate=False) else 0
M64 = collect(X36, G34, evaluate=False)[G34] if G34 in collect(X36, G34, evaluate=False) else 0
M65 = collect(X36, G35, evaluate=False)[G35] if G35 in collect(X36, G35, evaluate=False) else 0
M66 = collect(X36, G36, evaluate=False)[G36] if G36 in collect(X36, G36, evaluate=False) else 0

F_9_to_15 = Matrix([[M11, M12, M13, M14, M15, M16],
                    [M21, M22, M23, M24, M25, M26],
                    [M31, M32, M33, M34, M35, M36],
                    [M41, M42, M43, M44, M45, M46],
                    [M51, M52, M53, M54, M55, M56],
                    [M61, M62, M63, M64, M65, M66]])

B1 = X24.subs([(G24, 0), (G25, 0), (G26, 0), (G34, 0), (G35, 0), (G36, 0)])
B2 = X25.subs([(G24, 0), (G25, 0), (G26, 0), (G34, 0), (G35, 0), (G36, 0)])
B3 = X26.subs([(G24, 0), (G25, 0), (G26, 0), (G34, 0), (G35, 0), (G36, 0)])
B4 = X34.subs([(G24, 0), (G25, 0), (G26, 0), (G34, 0), (G35, 0), (G36, 0)])
B5 = X35.subs([(G24, 0), (G25, 0), (G26, 0), (G34, 0), (G35, 0), (G36, 0)])
B6 = X36.subs([(G24, 0), (G25, 0), (G26, 0), (G34, 0), (G35, 0), (G36, 0)])
B = - Matrix([[B1], [B2], [B3], [B4], [B5], [B6]])


print("\n=== It takes sympy too long to solve 6 linear equations analytically, but Maple does it well. ====")
print("\n=== Here we only solve numerically. ====")

F_9_to_15_num = F_9_to_15.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_), (G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num), (G22,G22_num), (G23,G23_num), (G33,G33_num)])

B_num = B.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_), (G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num), (G22,G22_num), (G23,G23_num), (G33,G33_num)])

G_mRNA_protein = F_9_to_15_num.inv()*B_num

G24_num, G25_num, G26_num, G34_num, G35_num, G36_num = G_mRNA_protein[0,0], G_mRNA_protein[1,0], G_mRNA_protein[2,0], G_mRNA_protein[3,0], G_mRNA_protein[4,0], G_mRNA_protein[5,0]

print("\n===== SOLVING EQUATIONS 16,17,18,19,20,21 FOR PROTEIN-PROTEIN CORRELATIONS AND VARIANCES ======")
# O2_16_17_18_19_20_21 = linsolve([X44, X45, X46, X55, X56, X66], G44, G45, G46, G55, G56, G66)


# G24 = O2_16_17_18_19_20_21.args[0][0]
# G25 = O2_16_17_18_19_20_21.args[0][1]
# G26 = O2_16_17_18_19_20_21.args[0][2]
# G34 = O2_16_17_18_19_20_21.args[0][3]
# G35 = O2_16_17_18_19_20_21.args[0][4]
# G36 = O2_16_17_18_19_20_21.args[0][5]

############# INVERTING MATRIX ############
M11 = collect(X44, G44, evaluate=False)[G44] if G44 in collect(X44, G44, evaluate=False) else 0
M12 = collect(X44, G45, evaluate=False)[G45] if G45 in collect(X44, G45, evaluate=False) else 0
M13 = collect(X44, G46, evaluate=False)[G46] if G46 in collect(X44, G46, evaluate=False) else 0
M14 = collect(X44, G55, evaluate=False)[G55] if G55 in collect(X44, G55, evaluate=False) else 0
M15 = collect(X44, G56, evaluate=False)[G56] if G56 in collect(X44, G56, evaluate=False) else 0
M16 = collect(X44, G66, evaluate=False)[G66] if G66 in collect(X44, G66, evaluate=False) else 0

M21 = collect(X45, G44, evaluate=False)[G44] if G44 in collect(X45, G44, evaluate=False) else 0
M22 = collect(X45, G45, evaluate=False)[G45] if G45 in collect(X45, G45, evaluate=False) else 0
M23 = collect(X45, G46, evaluate=False)[G46] if G46 in collect(X45, G46, evaluate=False) else 0
M24 = collect(X45, G55, evaluate=False)[G55] if G55 in collect(X45, G55, evaluate=False) else 0
M25 = collect(X45, G56, evaluate=False)[G56] if G56 in collect(X45, G56, evaluate=False) else 0
M26 = collect(X45, G66, evaluate=False)[G66] if G66 in collect(X45, G66, evaluate=False) else 0

M31 = collect(X46, G44, evaluate=False)[G44] if G44 in collect(X46, G44, evaluate=False) else 0
M32 = collect(X46, G45, evaluate=False)[G45] if G45 in collect(X46, G45, evaluate=False) else 0
M33 = collect(X46, G46, evaluate=False)[G46] if G46 in collect(X46, G46, evaluate=False) else 0
M34 = collect(X46, G55, evaluate=False)[G55] if G55 in collect(X46, G55, evaluate=False) else 0
M35 = collect(X46, G56, evaluate=False)[G56] if G56 in collect(X46, G56, evaluate=False) else 0
M36 = collect(X46, G66, evaluate=False)[G66] if G66 in collect(X46, G66, evaluate=False) else 0

M41 = collect(X55, G44, evaluate=False)[G44] if G44 in collect(X55, G44, evaluate=False) else 0
M42 = collect(X55, G45, evaluate=False)[G45] if G45 in collect(X55, G45, evaluate=False) else 0
M43 = collect(X55, G46, evaluate=False)[G46] if G46 in collect(X55, G46, evaluate=False) else 0
M44 = collect(X55, G55, evaluate=False)[G55] if G55 in collect(X55, G55, evaluate=False) else 0
M45 = collect(X55, G56, evaluate=False)[G56] if G56 in collect(X55, G56, evaluate=False) else 0
M46 = collect(X55, G66, evaluate=False)[G66] if G66 in collect(X55, G66, evaluate=False) else 0

M51 = collect(X56, G44, evaluate=False)[G44] if G44 in collect(X56, G44, evaluate=False) else 0
M52 = collect(X56, G45, evaluate=False)[G45] if G45 in collect(X56, G45, evaluate=False) else 0
M53 = collect(X56, G46, evaluate=False)[G46] if G46 in collect(X56, G46, evaluate=False) else 0
M54 = collect(X56, G55, evaluate=False)[G55] if G55 in collect(X56, G55, evaluate=False) else 0
M55 = collect(X56, G56, evaluate=False)[G56] if G56 in collect(X56, G56, evaluate=False) else 0
M56 = collect(X56, G66, evaluate=False)[G66] if G66 in collect(X56, G66, evaluate=False) else 0

M61 = collect(X66, G44, evaluate=False)[G44] if G44 in collect(X66, G44, evaluate=False) else 0
M62 = collect(X66, G45, evaluate=False)[G45] if G45 in collect(X66, G45, evaluate=False) else 0
M63 = collect(X66, G46, evaluate=False)[G46] if G46 in collect(X66, G46, evaluate=False) else 0
M64 = collect(X66, G55, evaluate=False)[G55] if G55 in collect(X66, G55, evaluate=False) else 0
M65 = collect(X66, G56, evaluate=False)[G56] if G56 in collect(X66, G56, evaluate=False) else 0
M66 = collect(X66, G66, evaluate=False)[G66] if G66 in collect(X66, G66, evaluate=False) else 0

F_16_to_21 = Matrix([[M11, M12, M13, M14, M15, M16],
                     [M21, M22, M23, M24, M25, M26],
                     [M31, M32, M33, M34, M35, M36],
                     [M41, M42, M43, M44, M45, M46],
                     [M51, M52, M53, M54, M55, M56],
                     [M61, M62, M63, M64, M65, M66]])

B1 = X44.subs([(G44, 0), (G45, 0), (G46, 0), (G55, 0), (G56, 0), (G66, 0)])
B2 = X45.subs([(G44, 0), (G45, 0), (G46, 0), (G55, 0), (G56, 0), (G66, 0)])
B3 = X46.subs([(G44, 0), (G45, 0), (G46, 0), (G55, 0), (G56, 0), (G66, 0)])
B4 = X55.subs([(G44, 0), (G45, 0), (G46, 0), (G55, 0), (G56, 0), (G66, 0)])
B5 = X56.subs([(G44, 0), (G45, 0), (G46, 0), (G55, 0), (G56, 0), (G66, 0)])
B6 = X66.subs([(G44, 0), (G45, 0), (G46, 0), (G55, 0), (G56, 0), (G66, 0)])
B = - Matrix([[B1], [B2], [B3], [B4], [B5], [B6]])

print("\n=== It takes sympy too long to solve 6 linear equations analytically, but Maple does it well. ====")
print("=== Here we only solve numerically. ====")
print("Note that due to the simplicity of the x6x6 equation the system of 6 equations\n can be reduced to 5+1 equations and solved analytically even with sympy.\n")

F_16_to_21_num = F_16_to_21.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_), (G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num), (G22,G22_num), (G23,G23_num), (G33,G33_num), (G24,G24_num), (G25,G25_num), (G26,G26_num), (G34,G34_num), (G35,G35_num), (G36,G36_num)])

B_num = B.subs([(N,N_), (l1p,l1p_), (l1m,l1m_), (l2,l2_), (theta01,theta01_), (theta10,theta10_), (eta11,eta11_), (gamma11,gamma11_), (k0,k0_), (k1,k1_), (nu01,nu01_), (nu10,nu10_), (tau1,tau1_), (tau2,tau2_), (tau3,tau3_), (Pon, Pon_), (G1,G1_num), (G2,G2_num), (G3,G3_num), (G4,G4_num), (G5,G5_num), (G6,G6_num), (G11, G11_num), (G12,G12_num), (G13,G13_num), (G14,G14_num), (G15,G15_num), (G16,G16_num), (G22,G22_num), (G23,G23_num), (G33,G33_num), (G24,G24_num), (G25,G25_num), (G26,G26_num), (G34,G34_num), (G35,G35_num), (G36,G36_num)])

G_protein_protein = F_16_to_21_num.inv()*B_num

G44_num, G45_num, G46_num, G55_num, G56_num, G66_num = G_protein_protein[0,0], G_protein_protein[1,0], G_protein_protein[2,0], G_protein_protein[3,0], G_protein_protein[4,0], G_protein_protein[5,0]

print("Protein variances:")
print("Soma: var(l_0)=<l_0^2> - <l_0>^2 = G44 - G4*(G4-1)", G44_num - G4_num*(G4_num-1))
print("Dendrite: var(l_1)=<l_1^2> - <l_1>^2 = G55 - G5*(G5-1)", G55_num - G5_num*(G5_num-1))
print("Spine: var(q_11)=<q_11^2> - <q_11>^2 = G66 - G6*(G6-1)", G66_num - G6_num*(G6_num-1))

print("\n**************** SUMMARY ****************\n")
print("\n========== MEANS: ==========")
print("<Active gene replicas>=", G1_num)
print("<mRNAs in soma>=", G2_num)
print("<mRNAs in dendrite>=", G3_num)
print("<Proteins in soma>=", G4_num)
print("<Proteins in dendrite>=", G5_num)
print("<Proteins in synapse>=", G6_num)

print("\n======= VARIANCES: ========")
print("RMS(Active gene replicas)=", sqrt(G11_num - G1_num*(G1_num-1)))
print("RMS(mRNAs in soma)=", sqrt(G22_num - G2_num*(G2_num-1)))
print("RMS(mRNAs in dendrite)=", sqrt(G33_num - G3_num*(G3_num-1)))
print("RMS(proteins in soma)=", sqrt(G44_num - G4_num*(G4_num-1)))
print("RMS(proteins in dendrite)=", sqrt(G55_num - G5_num*(G5_num-1)))
print("RMS(proteins in synapse)=", sqrt(G66_num - G6_num*(G6_num-1)))





# F_16_to_21 = Matrix([[-(theta01+1/tau3), theta10, 0, 0, 0],
#                      [theta01, -(eta11+theta01+theta10+2/tau3), gamma11, theta10, 0],
#                      [0, eta11, -(gamma11+theta01+1/tau3), 0, theta10],
#                      [0, theta01, 0, -(eta11+theta01+1/tau3), gamma11],
#                      [0, 0, theta01, eta11, -(gamma11+theta10+1/tau3)]])

# R = - Matrix([[k0*G24], [k0*G25+k1*G34], [k0*G26], [k1*G35], [k1*G36]])












































################# CHECKS #################
# G1_ = tau1*N*l1p
# G2_ = N*l1p*tau1*l2*tau2*(tau2*nu10+1)/(tau2*(nu01+nu10)+1)
# G3_ = N*l1p*tau1*l2*tau2**2*nu01/(tau2*(nu01+nu10)+1)
# # G4_ = (k0*(theta10*tau3+1)*(tau2*nu10+1)+k1*tau3*theta10*tau2*nu01)/((tau3*(theta01+theta10)+1)*(tau2*(nu01+nu10)+1))*N*l1p*tau1*l2*tau2*tau3
# G4_ = (k0*(theta10*tau3+1)*tau3*G2_ + k1*tau3**2*theta10*G3_)/((theta10+theta01)*tau3+1)
# # G5_ = (k1*nu01*tau2*(tau3*(theta01+2*theta10)+1) + theta01*k0*tau3*(theta10*tau3+1)*(tau2*nu10+1))/((tau3*(theta01+theta10)+1)*(tau2*(nu01+nu10)+1))*N*l1p*tau1*l2*tau2*tau3
# G5_ = tau3*k1/(theta10*tau3+1)*G3_ + theta01*tau3/(theta10*tau3+1)*G4_
# G6_ = eta11/gamma11*G5_

# print("-------------- CHECKS FIRST ORDER ------------------\n")
# print("_________________x1________________\n")
# pretty_print(cancel((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
# #print((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("_________________x2________________\n")
# pretty_print(cancel((collect(expand(RHS - LHS), x2, evaluate=False)[x2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
# print("\n")

# print("_________________x3________________\n")
# pretty_print(cancel((collect(expand(RHS - LHS), x3, evaluate=False)[x3]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
# print("\n")

# print("_________________x4________________\n")
# pretty_print(simplify((collect(expand(RHS - LHS), x4, evaluate=False)[x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
# print("\n")

# print("_________________x5________________\n")
# pretty_print(simplify((collect(expand(RHS - LHS), x5, evaluate=False)[x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
# print("\n")

# print("_________________x6________________\n")
# pretty_print(cancel((collect(expand(RHS - LHS), x6, evaluate=False)[x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
# print("\n")
