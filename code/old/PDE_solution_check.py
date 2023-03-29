import numpy as np
from sympy import *

# x1 = symbols('x1')
# t = symbols('t')
# N1 = symbols('N1')
# m = symbols('m')
# l1p = symbols('lambda1+')
# l1m = symbols('lambda1-')

# tau1 = (l1p+l1m)**-1
# Pon = l1p/(l1p+l1m)

# G = (x1*Pon*(1-exp(-t/tau1)) + 1)**N1*( x1*exp(-t/tau1)/(x1*Pon*(1-exp(-t/tau1))+1))**m

# RHS = l1p*N1*x1*G - (l1p*(x1+1)+l1m)*x1*diff(G, x1)

# print(simplify(RHS - diff(G,t)) == 0)


##### CHECKING VARIANCE FROM PAULSSON #####
N1 = symbols('N1')
n1 = symbols('n1')
n2 = symbols('n2')
n3 = symbols('n3')
l1p = symbols('lambda1^+')
l1m = symbols('lambda1^-')
l2 = symbols('lambda2')
l3 = symbols('lambda3')

tau1 = 1/(l1p+l1m)
tau2 = symbols('tau2')
tau3 = symbols('tau3')

Pon = l1p/(l1p+l1m)

G1 = Pon*N1
G2 = l2*tau2*Pon*N1
G3 = l3*tau3*l2*tau2*Pon*N1

G33 = tau2/(tau2+tau3)*(tau3*l3)**2*tau2*l2*Pon*N1* (1 + tau1*tau2*l2/(tau1+tau2)*(Pon*(N1-1)+l1p*N1*tau2+1) + l2*tau1*tau3/(tau1+tau3)*(l1p*N1*tau3+tau1/(tau1+tau2)*(N1*(l1p*tau2+Pon)+1-Pon)))

var3 = G33 + G3 - G3**2

var3_article = G3**2*(1/G3+1/G2*tau2/(tau2+tau3) + (1-Pon)/G1*tau2/(tau2+tau3)*tau1/(tau1+tau3)*(tau1+tau3+tau1*tau3/tau2)/(tau1+tau2))

pretty_print(simplify(var3/var3_article))
