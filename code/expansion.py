import numpy as np
from sympy import *

N = symbols('N')
n1 = symbols('n1')
n2 = symbols('n2')
n3 = symbols('n3')
l1p = symbols('lambda1^+')
l1m = symbols('lambda1^-')
l2 = symbols('lambda2')
l3 = symbols('lambda3')
theta01 = symbols('theta01')
theta10 = symbols('theta10')
eta11 = symbols('eta11')
gamma11 = symbols('gamma11')

k0 = symbols('kappa0')
k1 = symbols('kappa1')
nu01 = symbols('nu01')
nu10 = symbols('nu10')

tau1 = 1/(l1p+l1m)
tau2 = symbols('tau2')
tau3 = symbols('tau3')

Pon = l1p/(l1p+l1m)

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

# print("-------------- FIRST ORDER ------------------\n")
# print("_________________x1________________\n")
# pretty_print((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# #print((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("_________________x2________________\n")
# pretty_print((collect(expand(RHS - LHS), x2, evaluate=False)[x2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("_________________x3________________\n")
# pretty_print((collect(expand(RHS - LHS), x3, evaluate=False)[x3]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("_________________x4________________\n")
# pretty_print((collect(expand(RHS - LHS), x4, evaluate=False)[x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("_________________x5________________\n")
# pretty_print((collect(expand(RHS - LHS), x5, evaluate=False)[x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("_________________x6________________\n")
# pretty_print((collect(expand(RHS - LHS), x6, evaluate=False)[x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
# print("\n")

# print("-------------- SECOND ORDER ------------------\n")






################# CHECKS #################
G1_ = tau1*N*l1p
G2_ = N*l1p*tau1*l2*tau2*(tau2*nu10+1)/(tau2*(nu01+nu10)+1)
G3_ = N*l1p*tau1*l2*tau2**2*nu01/(tau2*(nu01+nu10)+1)
# G4_ = (k0*(theta10*tau3+1)*(tau2*nu10+1)+k1*tau3*theta10*tau2*nu01)/((tau3*(theta01+theta10)+1)*(tau2*(nu01+nu10)+1))*N*l1p*tau1*l2*tau2*tau3
G4_ = (k0*(theta10*tau3+1)*tau3*G2_ + k1*tau3**2*theta10*G3_)/((theta10+theta01)*tau3+1)
# G5_ = (k1*nu01*tau2*(tau3*(theta01+2*theta10)+1) + theta01*k0*tau3*(theta10*tau3+1)*(tau2*nu10+1))/((tau3*(theta01+theta10)+1)*(tau2*(nu01+nu10)+1))*N*l1p*tau1*l2*tau2*tau3
G5_ = tau3*k1/(theta10*tau3+1)*G3_ + theta01*tau3/(theta10*tau3+1)*G4_
G6_ = eta11/gamma11*G5_

print("-------------- CHECKS FIRST ORDER ------------------\n")
print("_________________x1________________\n")
pretty_print(cancel((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
#print((collect(expand(RHS - LHS), x1, evaluate=False)[x1]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0)]))
print("\n")

print("_________________x2________________\n")
pretty_print(cancel((collect(expand(RHS - LHS), x2, evaluate=False)[x2]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
print("\n")

print("_________________x3________________\n")
pretty_print(cancel((collect(expand(RHS - LHS), x3, evaluate=False)[x3]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
print("\n")

print("_________________x4________________\n")
pretty_print(simplify((collect(expand(RHS - LHS), x4, evaluate=False)[x4]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
print("\n")

print("_________________x5________________\n")
pretty_print(simplify((collect(expand(RHS - LHS), x5, evaluate=False)[x5]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
print("\n")

print("_________________x6________________\n")
pretty_print(cancel((collect(expand(RHS - LHS), x6, evaluate=False)[x6]).subs([(x1, 0), (x2, 0), (x3, 0), (x4, 0), (x5, 0), (x6, 0), (G1, G1_), (G2, G2_), (G3, G3_), (G4, G4_), (G5, G5_), (G6, G6_)])))
print("\n")
