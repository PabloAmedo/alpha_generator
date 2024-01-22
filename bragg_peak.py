# -*- coding: utf-8 -*-
"""
CreaEed on Wed Apr 19 11:08:56 2023

@auEhor: jacob
"""

import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd

rho_Ar = 1.784             #g/cm3
def bragg_peak(E,alpha_range):
  
    me=0.511                           #MeV/c**2
    M=931.5*4                          #MeV/c**2
    re=2.8179402894e-13                #cm
    NA=6.022e23                        #1/mol
    Co = 4*np.pi*NA*re**2*me           #MeV*cm2/mol (??)
    
    #-------Ar------
    I=(210)*1e-6                       #eV                                     #xq 210?? segun NIST 188
    A = 40                                                                     #Calculate Z/A effective factor for Ar/CF4
    Z = 18
    z = 2
    rho=1.66201e-3         #1.748*1000/(100**3)                                #according NIST rho=1.66201e-3 #g/cm3
    #---------------
      
    Sp=[]                           
    x=[]
    r=np.linspace(alpha_range/1000,alpha_range*10,int(1e5/2.5))
         
    for i in range(len(r)-1):
        gamma=E/M+1
        beta=(1-1/(gamma**2))**0.5
        dE_dx=rho*Co*Z/A*z**2/beta**2*(np.log(2*me/I)+np.log(beta**2/(1-beta**2))-beta**2)
        
        Sp.append(dE_dx)
        x.append((r[i+1]+r[i])/2)
        
        E=E-dE_dx*(r[i+1]-r[i])
        if x[i]>alpha_range and x[i]<alpha_range*1.1:
            break
    
    Sp=np.array(Sp)
    x=np.array(x)
    
    #Debemos cortar la distribución cuando baje a prácticamente cero
    #ya que hacemos random*Rcsda 
    
    
    # plt.figure(figsize=(8,7),dpi=100)
    # plt.grid(True)
    # plt.ylabel('Stopping power/Alpha energy (1/cm)')
    # plt.xlabel('x(cm)')
    # # plt.xlim(0,6)
    # plt.plot(x,Sp/np.sum(Sp))
    # print('x',x)
    # print('Sp',Sp)
    
    x=x/x[-1] # 

    # y=Sp/np.sum(Sp)
    
    acum=[]
    acum=np.cumsum(Sp)/np.sum(Sp)
    acum=np.array(acum)    # Integral // valor max=1
    
    return (x,Sp,acum)

def random_bragg(x,acum):
    '''
    It gives a random number following a given distribution
    '''
    aleat=np.random.rand() # Numero aleatorio 0 a 1
    
    diff=acum-aleat
    diff=np.absolute(diff)
    min_diff=np.min(diff)
        
    i=np.where(diff==min_diff)
    I=i[0]

    return x[I]

"""
# TEST

# x,Sp,acum=bragg_peak(5.5,5)
# y=0
# i=0
# mayor=[]
# menor=[]
# values=[]
# plt.figure()
# for i in range(10000):
#     i=i+1
#     y=random_bragg(x,acum)
#     if y>0.7:
#         color='red'
#         mayor.append(y)
#     else:
#         color='blue'
#         menor.append(y)
#     plt.plot(i,y,'.',color=color) 
#     values.append(y)
    
# hist, bins = np.histogram(values,bins=100)
# plt.figure()
# plt.bar(bins[:-1], hist, width=np.diff(bins)*0.9, align='edge')
# plt.show()



# Dio un 50% los puntos mayores a 0.8
"""
NIST = pd.read_csv('C:/Users/usuario/Desktop/NIST-ASTAR.txt' ,delimiter='\t')
x, Sp, acum = bragg_peak(5.5, 4.5)
x_N = np.linspace(min(x), max(x), len(NIST))
plt.figure()
plt.title('Bragg peak', fontsize = 20)
plt.plot(x*4.5, Sp*rho_Ar, '.', color = 'black', label = 'simulated points')

plt.xlabel('length [cm]')
plt.ylabel('stopping power [MeV $cm^{-1}$]')
plt.grid()
plt.legend()
