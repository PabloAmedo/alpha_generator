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
gammas_Ar=np.array([4.0,3.5])
cluster_Ar=np.array([27.8,28.6])

gammas_Xe = np.array(( 4.0))
clusters_Xe = np.array((44))

Emuon = 4000
def clusters_cm(Emuon, Wi = 26.4, gas = 'Argon', cl_measurements = cluster_Ar, gammas = gammas_Ar):
    
    folder= 'NIST data/'
    
    NISTdata=pd.read_csv(folder + gas + '_Sp.txt', delimiter='\t', header = 4)

    NISTdata.rename(columns={'MeV      ':'Energy [MeV]',
                             'MeV cm2/g':'dE/dx [MeV cm2/g]'},
                    inplace=True)

    NISTdata=NISTdata.drop(columns={'Unnamed: 2'})

    #Manually introduce experimental data (from book)  --> This is just to make the W calculation (we need data)
    
    E_book = 0.511*(gammas - 1)                                                    #transform gamma to E
    
    
    
    #K calculation----------------------------------------------------------------- 
    
    """
    We are going to calculate the difference between the data points ([ERM69]) and 
    the corresponding stopping power predicted by Bethe-Bloch curve:
    
    In a first approach I am going to perform a linear extrapolation for the BB 
    curve to match the energies, but I can consider two more options for a future:
        FIXME
        
        -Export NIST data with the relevant energies included (they do the extrapo)
        
        -Fit a function to NIST data. (I guess that this is more general but I am 
                                       too lazy for it)
    """
    #Localize the energy range where we are moving on for every data point
    
    #Linear fit through it
    def linear(x,a,b):
        return a*x+b
    
    if gas == 'Argon':
        E1_a=NISTdata.loc[NISTdata['Energy [MeV]']<E_book[0]].index[-1]
        E1_b=(NISTdata.loc[NISTdata['Energy [MeV]']>E_book[0]]).index[0]
        
        E2_a=NISTdata.loc[NISTdata['Energy [MeV]']<E_book[1]].index[-1]
        E2_b=(NISTdata.loc[NISTdata['Energy [MeV]']>E_book[1]]).index[0]
    
            
        
        
        xdata=[NISTdata.loc[E1_a]['Energy [MeV]'], NISTdata.loc[E1_b]['Energy [MeV]']]#, NISTdata.loc[E2_a]['Energy [MeV]']]
        ydata=[NISTdata.loc[E1_a]['dE/dx [MeV cm2/g]'], NISTdata.loc[E1_b]['dE/dx [MeV cm2/g]']]#, NISTdata.loc[E2_a]['dE/dx [MeV cm2/g]']]
        fit,_=sco.curve_fit(linear, xdata, ydata)
        
        #Calculate dE/dx for data points energy
        dEdx_fit=linear(E_book,*fit)
        
        K = cl_measurements / dEdx_fit
        Kmean=np.mean(K)
        
    elif gas == 'Xenon':
        
        E1_a=NISTdata.loc[NISTdata['Energy [MeV]']<E_book].index[-1]
        E1_b=NISTdata.loc[NISTdata['Energy [MeV]']>E_book].index[0]
        
        xdata=[NISTdata.loc[E1_a]['Energy [MeV]'], NISTdata.loc[E1_b]['Energy [MeV]']] #, NISTdata.loc[E2_a]['Energy [MeV]']]
        ydata=[NISTdata.loc[E1_a]['dE/dx [MeV cm2/g]'], NISTdata.loc[E1_b]['dE/dx [MeV cm2/g]']]#, NISTdata.loc[E2_a]['dE/dx [MeV cm2/g]']]
        fit,_=sco.curve_fit(linear, xdata, ydata)
        
        #Calculate dE/dx for data points energy
        dEdx_fit=linear(E_book,*fit)
        
        K = cl_measurements / dEdx_fit
        Kmean=np.mean(K)
        
        
        
        
    """
    print('\n\n')
    print('RESULTS for the K factor in {}:'.format(gas))
    print('* K_i factor:\t', K)
    print('* K mean:\t\t', Kmean)
    print('-'*59)
    """
    #------------------------------------------------------------------------------
    #CALCULATION OF THE AVG N OF CLUSTERS ACCORDING TO MUON ENERGY
    
    #This energy is higher than the data used from NIST, so I am just extrapolating
    #extrapolate to linear y=mx+n
    m_extr=(NISTdata['dE/dx [MeV cm2/g]'][80]-NISTdata['dE/dx [MeV cm2/g]'][78])/(NISTdata['Energy [MeV]'][80]-NISTdata['Energy [MeV]'][78])
    n_extr=NISTdata['dE/dx [MeV cm2/g]'][80]-m_extr*NISTdata['Energy [MeV]'][80]
    
    Ncluster_muon = (m_extr*Emuon+n_extr)*Kmean
    print('\nFor a muon with {} MeV:\t'.format(Emuon), Ncluster_muon,'[clusters/cm]')
    
    return (m_extr*Emuon+n_extr) , Ncluster_muon

#PLOT--------------------------------------------------------------------------
"""
#THIS CAN BE AVOIDED IN THE FUNCTION -- I USED IT FOR SOME MASTER'S TASKS

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

