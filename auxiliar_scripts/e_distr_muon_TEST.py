# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 18:11:46 2023

@author: diego
"""

import numpy as np

print('Running...')
#<n_{e, cl}> CALCULATION

W_i         =   26.4                                                           #eV
rho_Ar      =   1.662e-3                                                       #g/cm3
#next quantities taken from Clusters_calculation.py 
dEdx_4GeV   =   3.25                                                            #MeV cm2/g
Ncl         =   29.47 * 1                                                           #clusters/cm

probs_less_20=np.array([65.6,15,6.4,3.5,2.25,1.55,1.05,0.81,0.61,0.49,0.39,0.30,
   0.25,0.2,0.16,0.12,0.095,0.075,0.063])/100

n_prom = ( ( dEdx_4GeV ) * rho_Ar ) / ( Ncl * W_i * (1e-6) )

#------------------------------------------------------------------------------

#n_max CALCULATION --> This has to fit the  n_max

n = 1e8
num = 0
denom = 0

for i in range(1,int(n)):
    if i <= 19:
        num = num + probs_less_20[i-1]
        denom = denom + ( probs_less_20[i-1] / i**2 )
    else:
        num = num + 0.216 / i
        denom = denom + ( 0.216 / i**2 )
    
    if num/denom > n_prom:
        print('n_max\t=\t', i, '->\t', num/denom)
        break
    
n_max = num/denom
#------------------------------------------------------------------------------

#MAX E TO DELTA e-

gamma = 4
beta =  0.9996
m_e = 0.511                     #MeV/c2
M_muon = 105.66                 #MeV/c2

W_max = (2 * m_e * (beta * gamma)**2)

print('W_max\t\t=\t', W_max )

