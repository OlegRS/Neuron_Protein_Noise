PC = np.zeros([15,15]) #Pearson correlation coefficient
for i in range(0, 15):
    PC[i,i] = 1
    for j in range(i+1, 15):
        PC[i,j] = (O2[conv[i,i]] - E_num[i]*E_num[j])/(RMS[i]*RMS[j])
        PC[j,i] = PC[i,j]
print('Pearson correlation matrix:\n', PC)


# PC = np.zeros([15,15]) #Pearson correlation coefficient
# for i in range(0, 15):
#     for j in range(i+1, 15):
#         PC[i,j] = (O2[conv[i,j]] - E[i]*E[j])/(RMSs[i]*RMSs[j])

# RHS = np.array(RHS, dtype=float)
# XXX = np.array(XXX, dtype=float)
# print('Inverting the matrix...')
# XXX_inv = np.linalg.inv(XXX)
# O2 = np.dot(XXX_inv, RHS)

# RMSs = np.zeros(15)
# # print('Second order:\n', O2)

# print('RMSs:\n')
# E_num = np.array(E_num)
# for i in range(0,15):
#     RMSs[i] = sqrt(O2[conv[i,i]] - (E_num[i] - 1)*E_num[i])

# SN_ratios = [sqrt(O2[conv[i,i]] - (E_num[i] - 1)*E_num[i])/E_num[i]*100 for i in range(0,15)] # Signal to noise ratios %

# print('Signal to Noise ratios %:', SN_ratios)


############################################################################################
# k = 0
# RHS = Matrix(symbols('X0:120'))
# conv = np.zeros([15, 15], dtype=int)
# for i_ in range(0,15):
#     for j_ in range(i_,15):
#         RHS[k] = -XX[k].subs([(x[j],0) for j in range(0,15)]
#                            +[(GG,0)]
#                            +[(GM[j],0) for j in range(0,len(GM))]
#                            +[(GP[j],0) for j in range(0,len(GP))]
#                            +[(MM[j],0) for j in range(0,len(MM))]
#                            +[(MP[j],0) for j in range(0,len(MP))]
#                            +[(PP[j],0) for j in range(0,len(PP))]
#                            )
#         conv[i_,j_]=k
#         conv[j_,i_]=k        
#         k+=1        
        
# RHS = np.array(RHS, dtype=float)
# XXX_ = np.array(XXX_, dtype=float)
# print('Inverting the matrix...')
# XXX_inv = np.linalg.inv(XXX_)
# O2 = np.dot(XXX_inv, RHS)

# RMSs = np.zeros(15)
# # print('Second order:\n', O2)

# print('RMSs:\n')
# E_num = np.array(E_num)
# for i in range(0,15):
#     RMSs[i] = sqrt(O2[conv[i,i]] - (E_num[i] - 1)*E_num[i])

# SN_ratios = [sqrt(O2[conv[i,i]] - (E_num[i] - 1)*E_num[i])/E_num[i]*100 for i in range(0,15)] # Signal to noise ratios %

############################################################################################





# for i in range(0,15):
#     for j in range(0,15):
#         Y[i] = collect(X[i], E[j])



# E[9] =  E[6]*eta11/gamma11
# E[10] = E[6]*eta12/gamma12
# E[11] = E[7]*eta21/gamma21
# E[12] = E[7]*eta22/gamma22
# E[13] = E[8]*eta31/gamma31
# E[14] = E[8]*eta32/gamma32


########## CREATING A MATRIX ###########



# Matrix([[-E0*lambda1^+ - E0*lambda1^- + N*lambda1^+],
#         [E0*lambda2 - E1*nu01 - E1/tau2 + E2*nu10],
#         [E1*nu01 - E2*nu10 - E2*nu12 - E2*nu13 - E2/tau2 + E3*nu21 + E4*nu31],
#         [E2*nu12 - E3*nu21 - E3/tau2],
#         [E2*nu13 - E4*nu31 - E4/tau2],
#         [E1*kappa0 - E5*theta01 - E5/tau3 + E6*theta10],
#         [E10*gamma12 + E2*kappa1 + E5*theta01 - E6*eta11 - E6*eta12 - E6*theta10 - E6*theta12 - E6*theta13 - E6/tau3 + E7*theta21 + E8*theta31 + E9*gamma11],
#         [E11*gamma21 + E12*gamma22 + E3*kappa2 + E6*theta12 - E7*eta21 - E7*eta22 - E7*theta21 - E7/tau3],
#         [E13*gamma31 + E14*gamma32 + E4*kappa3 + E6*theta13 - E8*eta31 - E8*eta32 - E8*theta31],
#         [E6*eta11 - E9*gamma11],
#         [-E10*gamma12 + E6*eta12],
#         [-E11*gamma21 + E7*eta21],
#         [-E12*gamma22 + E7*eta22],
#         [-E13*gamma31 + E8*eta31],
#         [-E14*gamma32 + E8*eta32]])




# Matrix([[-E0*lambda1^+ - E0*lambda1^- + N*lambda1^+],
#         [E0*lambda2 - E1*nu01 - E1/tau2 + E2*nu10],
#         [E1*nu01 - E2*nu10 - E2*nu12 - E2*nu13 - E2/tau2 + E3*nu21 + E4*nu31],
#         [E2*nu12 - E3*nu21 - E3/tau2],
#         [E2*nu13 - E4*nu31 - E4/tau2],
#         [E1*kappa0 - E5*theta01 - E5/tau3 + E6*theta10],
#         [E10*gamma12 + E2*kappa1 + E5*theta01 - E6*eta11 - E6*eta12 - E6*theta10 - E6*theta12 - E6*theta13 - E6/tau3 + E7*theta21 + E8*theta31 + E9*gamma11],
#         [E11*gamma21 + E12*gamma22 + E3*kappa2 + E6*theta12 - E7*eta21 - E7*eta22 - E7*theta21 - E7/tau3],
#         [E13*gamma31 + E14*gamma32 + E4*kappa3 + E6*theta13 - E8*eta31 - E8*eta32 - E8*theta31],
#         [E6*eta11 - E9*gamma11],
#         [-E10*gamma12 + E6*eta12],
#         [-E11*gamma21 + E7*eta21],
#         [-E12*gamma22 + E7*eta22],
#         [-E13*gamma31 + E8*eta31],
#         [-E14*gamma32 + E8*eta32]])
