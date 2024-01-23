# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:42:22 2023

@author: diego
"""

import numpy as np
import pandas as pd
import random as rand
import matplotlib.pyplot as plt
import scipy.optimize as sco

# CALCULATION OF CLUSTERS/MM 
plt.close('all')

"""
Notes about notation:
    
book: reference to the BlumRieglerRolandi text (sometimes I will directly reffer to its 
                                     bibliography as [] )
NIST: ESTAR data from NIST (well, ESTAR in this case, I can try to use ASTAR 2)
"""

#Load data and format----------------------------------------------------------
folder= 'NIST data/'
gas='Argon'         #gas name to choose the _NIST data
NISTdata=pd.read_csv(folder + gas + '_Sp.txt', delimiter='\t', header=4)

NISTdata.rename(columns={'MeV      ':'Energy [MeV]',
                         'MeV cm2/g':'dE/dx [MeV cm2/g]'},
                inplace=True)

NISTdata=NISTdata.drop(columns={'Unnamed: 2'})

#Manually introduce experimental data (from book)  --> This is just to make the W calculation (we need data)
gamma_book=np.array([4.0,3.5])
E_book=0.511*(gamma_book-1)                                                    #transform gamma to E
Ncluster=np.array([27.8,28.6])

Emuon = 500 #MeV
We = 25.7e-6 #Ionization energy for Ar

def clusters_cm(Emuon = Emuon, Wi = We, gas = 'Argon', exp_data = 'Argon'):
    #K calculation----------------------------------------------------------------- 
    
    """
    We are going to calculate the difference between the data points ([ERM69]) and 
    the corresponding stopping power predicted by Bethe-Bloch curve:
    
    In a first approach I am going to perform a linear extrapolation for the BB 
    curve to match the energies, but I can consider two more options for a future:
        *FIXME*
        
        -Export NIST data with the relevant energies included (they do the extrapo)
        
        -Fit a function to NIST data. (I guess that this is more general but I am 
                                       too lazy for it)
    """
    #Localize the energy range where we are moving on for every data point
    E1_a=NISTdata.loc[NISTdata['Energy [MeV]']<E_book[0]].index[-1]
    E1_b=(NISTdata.loc[NISTdata['Energy [MeV]']>E_book[0]]).index[0]
    
    E2_a=NISTdata.loc[NISTdata['Energy [MeV]']<E_book[1]].index[-1]
    E2_b=(NISTdata.loc[NISTdata['Energy [MeV]']>E_book[1]]).index[0]
    
    #Linear fit through it
    def linear(x,a,b):
        return a*x+b
    
    xdata=[NISTdata.loc[E1_a]['Energy [MeV]'], NISTdata.loc[E1_b]['Energy [MeV]']]#, NISTdata.loc[E2_a]['Energy [MeV]']]
    ydata=[NISTdata.loc[E1_a]['dE/dx [MeV cm2/g]'], NISTdata.loc[E1_b]['dE/dx [MeV cm2/g]']]#, NISTdata.loc[E2_a]['dE/dx [MeV cm2/g]']]
    fit,_=sco.curve_fit(linear, xdata, ydata)
    
    #Calculate dE/dx for data points energy
    dEdx_fit=linear(E_book,*fit)
    
    K=Ncluster/dEdx_fit
    Kmean=np.mean(K)
    print('\n\n')
    print('RESULTS for the K factor in {}:'.format(gas))
    print('* K_i factor:\t', K)
    print('* K mean:\t\t', Kmean)
    print('-'*59)
    
    
    
    
    #------------------------------------------------------------------------------
    #CALCULATION OF THE AVG N OF CLUSTERS ACCORDING TO MUON ENERGY
    
    """
    FIXME: I am setting a fix energy to the muon just to test. This has to be 
    updated to a E distribution.
    
    I can make here a try-except where if the E is in the range of NIST data we do 
    the same as before, and if it is not we can extrapolate as a linear regression. 
    I mean, it won't be too much larger than a few GeV... it should work
    """
    
    #This energy is higher than the data used from NIST, so I am just extrapolating
    #extrapolate to linear y=mx+n
    m_extr=(NISTdata['dE/dx [MeV cm2/g]'][80]-NISTdata['dE/dx [MeV cm2/g]'][78])/(NISTdata['Energy [MeV]'][80]-NISTdata['Energy [MeV]'][78])
    n_extr=NISTdata['dE/dx [MeV cm2/g]'][80]-m_extr*NISTdata['Energy [MeV]'][80]
    
    Ncluster_muon = (m_extr*Emuon+n_extr)*Kmean
    print('\nFor a muon with {} MeV:\t'.format(Emuon), Ncluster_muon,'[clusters/mm]')
    
    return (m_extr*Emuon+n_extr) , Ncluster_muon
#------------------------------------------------------------------------------
#NUMBER ELECTRON PER CLUSTER AND PROBABILITY
"""
#THIS CAN BE AVOIDED IN THE FUNCTION -- IT'S DONE BY ANOTHER CLASS (?)

k=np.linspace(1,32,32) #n of e- (from BlumRieglerRolandi book)[FIS91]
probAr_data=np.array([65.6,15.0,6.4,3.5,2.25,1.55,1.05,0.81,0.61,0.49,0.39,0.3,
                 0.25,0.2,0.16,0.12,0.095,0.075,0.063]) #from here on it goes as ~1/n**2
probAr_extr=[]
for i in range(len(probAr_data), len(k)):
    probAr_extr.append(21.6/k[i]**2)

probAr_extr=np.array(probAr_extr)
probAr=np.concatenate((probAr_data,probAr_extr))
probAr=np.append(probAr,(100-sum(probAr)))
k=np.append(k,(len(k)+1))

#LA ULTIMA PROBABILIDAD ESTA PUESTA POR NORMALIZACION
acum=[]
#e- per cluster choice // random weighted
for i in range(10000):
    n_e=np.random.choice(a=k,p=probAr/100)
    acum.append(n_e)
"""

#PLOT--------------------------------------------------------------------------
"""
#THIS CAN BE AVOIDED IN THE FUNCTION -- I DID IT TO SHOW SOME FIGURES

fig, ax1=plt.subplots()

ax1.plot(NISTdata['Energy [MeV]'], NISTdata['dE/dx [MeV cm2/g]'], color='b', 
         label='NIST data')
ax1.set_xscale('log')


ax1.set_xlabel('E [MeV]', fontsize = 13)
ax1.set_ylabel('-dE/dx [MeVg$^{-1}cm^2$]', color='b', fontsize = 13)
ax1.set_title('Ar-gas', fontsize = 20)

ax2=ax1.twinx()

ax2.errorbar(E_book, Ncluster, yerr=[0.3,0.5], fmt='.',color='black', 
             label='[ERM69] data')
ax2.set_ylabel('N cluster/mm for {}'.format(gas), color='g', fontsize = 13)
ax2.set_ylim([0,max(Ncluster)+Kmean])
ax2.set_xscale('log')

ax2.plot(NISTdata['Energy [MeV]'], Kmean*NISTdata['dE/dx [MeV cm2/g]'],
         '--', color='g', label='NIST x{:.2f}'.format(Kmean))
ax2.plot(Emuon, Ncluster_muon, '*', color='r', markersize = 10, label = '$E_{\mu}$')

ax2.spines['left'].set_color('b')
ax2.spines['right'].set_color('g')
fig.legend(loc = 'lower left')
fig.tight_layout()
"""

