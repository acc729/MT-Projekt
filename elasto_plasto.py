# -- coding: utf-8 --
"""
Created on Mon Aug 17 19:08:27 2020

@author: Lennard
"""


import numpy as np
# import main
import derf2
import tensors






#Material Parameter (0=E; 1=mue; 2=H; 3=M; 4=a; 5=b; 6=h; 7=G)
mat_pams =np.array([69000000,0.33,2000000,8,1.24,1.02,1.15,25.510**6])
delta = 1e-6
tol = 1e-6
res = tol+1

kappa_n = 0 # has to be input
I = np.eye(3)

#stress state 0=sig11; 1=sig22; 2=sig12


   
eps = np.zeros(6)      #has to be input 

der_f = np.zeros(6)
N = np.zeros(6)
A = np.zeros((7,7))
R = np.zeros((7))
R_ba = np.zeros(8)
der_N = np.zeros((6,6)) 

#-------------------------------------
ivar_n = np.zeros(6)     #
#---------------

eps = np.array([1, 0, 0, 0, 0, 0])      #input!

# sig = np.array([10000, 0, 2000000, 0, 0, 0])



#yield stresses and function
sig_yield11=50*10**6
sig_yield22=sig_yield11/mat_pams[6]

i = 0



#tangent modulus

IdyI, Isym = tensors.identityTensors()
IsymDev = Isym - 1/3 * (IdyI)                          #ISYM ISYMDEV NOCHMAL PRÜFEN BITTE!!!!! 

G = (1/(2*(1+mat_pams[1])))*mat_pams[0]

K = (mat_pams[0]/(3-6*mat_pams[1]))

E_t = tensors.kelvin(2*G*IsymDev + K * IdyI)
sig_trial = np.inner(E_t, eps-ivar_n)


#trial function
K1 = (sig_trial[0] + mat_pams[6]*sig_trial[1])/2
K2 = np.sqrt(((sig_trial[0]-mat_pams[6]*sig_trial[1])/2)**2 + (mat_pams[5]*2)*(sig_trial[5])**2)

f = mat_pams[4]*(abs(K1+K2))**mat_pams[3] + mat_pams[4]*(abs(K1-K2))**mat_pams[3] + (2-mat_pams[4])*(abs(2*K2))**mat_pams[3]
g = (1/2*f)**(1/mat_pams[3])


phi = g-sig_yield11 #yield function

if phi < 0: #elastic case
    sig = sig_trial
    ivar_np1 = ivar_n 

else: #plastic case
    dgamma = 0
    #sig_np1 = np.zeros(6)
    kappa_n = sig_yield11
    kappa_np1 = kappa_n
    
    sig = 1*sig_trial
    

    while res > tol:
        
        K1 = (sig[0] + mat_pams[6]*sig[1])/2
        K2 = np.sqrt(((sig[0]-mat_pams[6]*sig[1])/2)**2 + (mat_pams[5]*2)*(sig[5])**2)
        
        
        f = mat_pams[4]*(abs(K1+K2))**mat_pams[3] + mat_pams[4]*(abs(K1-K2))**mat_pams[3] + (2-mat_pams[4])*(abs(2*K2))**mat_pams[3]
        g = (1/2*f)**(1/mat_pams[3])
        
        phi = g-kappa_np1             #yield function
        
        
        der_f[0] = mat_pams[3]/2*((mat_pams[4]*(K1-K2)*np.abs(K1-K2)**(mat_pams[3]-2)*(1-((sig[0]-mat_pams[6]*sig[1])/(2*K2)))+mat_pams[4]*(K1+K2)*np.abs(K1+K2)**(mat_pams[3]-2)*(1+((sig[0]-mat_pams[6]*sig[1])/(2*K2)))+
        2**mat_pams[3]*(2-mat_pams[4])*K2**(mat_pams[3]-1)*(sig[0]-mat_pams[6]*sig[1])/(2*K2)))
        
        der_f[1] = mat_pams[3]*mat_pams[6]/2*((mat_pams[4]*(K1+K2))*np.abs(K1+K2)**(mat_pams[3]-2)*(1-((sig[0]-mat_pams[6]*sig[1])/(2*K2)))+mat_pams[4]*(K1-K2)*np.abs(K1-K2)**(mat_pams[3]-2)*(1+((sig[0]-mat_pams[6]*sig[1])/(2*K2)))-
        2**mat_pams[3]*(2-mat_pams[4])*K2**(mat_pams[3]-1)*(sig[0]-mat_pams[6]*sig[1])/(2*K2))   
    
        der_f[5] = mat_pams[3]*mat_pams[5]**2*sig[5]/(1*K2)*(mat_pams[4]*(K1+K2)*np.abs(K1+K2)**(mat_pams[3]-2)-mat_pams[4]*(K1-K2)*np.abs(K1-K2)**(mat_pams[3]-2) + 2*(2*mat_pams[4]*(2*K2)**(mat_pams[3]-1)))
        
        
        
        N = 1/(2*mat_pams[3])*(f/2)**((1-mat_pams[3])/mat_pams[3])*der_f 
        
        
        # second derivatives of f wrt sigma
        ddf = derf2.der(mat_pams, sig)    
        
        
        # derivatives of N wrt sigma
        der_N =  1/(2*mat_pams[3])*(f/2)**((1-mat_pams[3])/mat_pams[3])*ddf+ (1-mat_pams[3])/(4*mat_pams[3]**2)*(f/2)**((1-2*mat_pams[3])/mat_pams[3])*(np.outer(der_f,der_f))
        
        
        A[0:6, 0:6] = np.linalg.inv((np.linalg.inv(E_t))+dgamma*der_N)
        # A[6,0] = 0
        # A[0,6] = 0
        A[6,6] = mat_pams[2]
        
        
        # residuum matrix as R
        
        R[0:6] = np.linalg.inv(E_t).dot((sig-sig_trial))+dgamma*N
        R[6] = 1/mat_pams[2]*(kappa_np1-kappa_n) - dgamma
        
        R_ba[0:7] = R
        R_ba[7] = phi
        
        C = np.zeros(7)
        C[0:6] = N
        C[6] = -1
        
        ddgamma = (phi - np.inner(C,np.inner(A,R)))/np.inner(np.inner(C,A), np.transpose(C))
        
        R1 = np.inner(-A, (R + ddgamma*np.transpose(C)))
        
        dgamma = dgamma + ddgamma
        sig = R1[0:6] + sig              
        kappa_np1 = R1[6] + kappa_np1
        
        res = np.linalg.norm(R_ba)
        
        
        i = i+1



    
                

                                                                                                                                                   