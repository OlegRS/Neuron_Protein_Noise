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


# print("============ SOLVING FIRST ORDER... ===========")

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


print("============ SOLVING EQUATIONS 1,2,3 ==============")
O2_1_2_3 = linsolve([X11, X12, X13], G11, G12, G13)

G11 = O2_1_2_3.args[0][0]
G12 = O2_1_2_3.args[0][1]
G13 = O2_1_2_3.args[0][2]

print("============ SOLVING EQUATIONS 4,5,6 ==============")
O2_4_5_6 = linsolve([X14, X15, X16], G14, G15, G16)

G14 = O2_4_5_6.args[0][0]
G15 = O2_4_5_6.args[0][1]
G16 = O2_4_5_6.args[0][2]


print("============ SOLVING EQUATIONS 7,8,12 ==============")
O2_7_8_12 = linsolve([X22, X23, X33], G22, G23, G33)

G22 = O2_7_8_12.args[0][0]
G23 = O2_7_8_12.args[0][1]
G33 = O2_7_8_12.args[0][2]

print("============ SOLVING EQUATIONS 9,10,11,13,14,15 ==============")
O2_9_10_11_13_14_15 = linsolve([X24, X25, X26, X34, X35, X36], G24, G25, G26, G34, G35, G36)

G24 = O2_9_10_11_13_14_15.args[0][0]
G25 = O2_9_10_11_13_14_15.args[0][1]
G26 = O2_9_10_11_13_14_15.args[0][2]
G34 = O2_9_10_11_13_14_15.args[0][3]
G35 = O2_9_10_11_13_14_15.args[0][4]
G36 = O2_9_10_11_13_14_15.args[0][5]


print("============ SOLVING EQUATIONS 9,10,11,13,14,15 ==============")
O2_16_17_18_19_20_21 = linsolve([X44, X45, X46, X55, X56, X66], G44, G45, G46, G55, G56, G66)

G24 = O2_16_17_18_19_20_21.args[0][0]
G25 = O2_16_17_18_19_20_21.args[0][1]
G26 = O2_16_17_18_19_20_21.args[0][2]
G34 = O2_16_17_18_19_20_21.args[0][3]
G35 = O2_16_17_18_19_20_21.args[0][4]
G36 = O2_16_17_18_19_20_21.args[0][5]


# print("============ SOLVING NEXT 17 EQUATIONS... ==============")

# O2_next_17 = linsolve([X14, X15, X16, X22, X23, X24, X25, X26, X33, X34, X35, X36, X44, X45, X46, X55, X56],
#                       G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56)
                      
# O2 = linsolve([X11, X12, X13, X14, X15, X16, X22, X23, X24, X25, X26, X33, X34, X35, X36, X44, X45, X46, X55, X56, X66],
#       G11, G12, G13, G14, G15, G16, G22, G23, G24, G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66)



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