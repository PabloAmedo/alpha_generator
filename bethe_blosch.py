# -*- coding: uEf-8 -*-
"""
CreaEed on Wed Apr 19 11:08:56 2023

@auEhor: jacob
"""
from math import*
import matplotlib.pyplot as plt
import numpy as np
import random

def bragg_peak(E,alpha_range):
  
    me=0.511     #MeV/c**2
    M=933        #MeV/c**2
    re=2.8e-13   #cm
    NA=6.02e23   #1/mol
    Co = 4*np.pi*NA*re**2*me
    
    #-------Ar------
    I=(15)*1e-6
    A = 40
    Z = 18
    z = 2
    rho=1.748*1000/(100**3)
    #---------------
      
    Sp=[]
    x=[]
    r=np.linspace(alpha_range/1000,alpha_range*10,int(1e5/2.5))
     
    for i in range(len(r)-1):
        gamma=E/M+1
        beta=(1-1/(gamma**2))**0.5
        dE_dx=rho*Co*Z/A*z**2/beta**2*(np.log(2*me/I)+np.log(beta**2/(1-beta**2)))  
        
        Sp.append(dE_dx)
        x.append((r[i+1]+r[i])/2)
        
        E=E-dE_dx*(r[i+1]-r[i])
        if E<0:
            break
    
    Sp=np.array(Sp)
    x=np.array(x)
    x=x/x[-1] # 
    
    # plt.figure(figsize=(8,7),dpi=100)
    # plt.grid(True)
    # plt.ylabel('Energy loss / Alphas energy')
    # plt.plot(x,Sp/np.sum(Sp))
    # print('x',x)
    # print('Sp',Sp)
    return x,Sp

def random_bragg(x,Sp):
    
    acum=[]
    for i in range(len(Sp)):
        acum.append(np.sum(Sp[0:i])/np.sum(Sp))
    # plt.figure()
    # plt.plot(x,acum)
    
    acum=np.array(acum)    # Integral // valor max=1
    aleat=np.random.rand() # Numero aleatorio 0 a 1
    
    diff=acum-aleat
    diff=np.absolute(diff)
    min_diff=np.min(diff)
        
    i=np.where(diff==min_diff)
    I=i[0]
    # print(i)
    # print(I)
    # print(x) #da 13????

    return x[I]


# TEST

# acum,diff,Sp,y=bragg_peak(5.5,5)
# y=0
# i=0
# x=[]
# mayora08=[]
# menor=[]
# for i in range(100):
#     i=i+1
#     acum,diff,Sp,y=bragg_peak(5.5,5)
#     if y>0.8:
#         color='red'
#         mayora08.append(y)
#     else:
#         color='blue'
#         menor.append(y)
#     plt.plot(i,y,'.',color=color)

# Dio un 50% los puntos mayores a 0.8
    

