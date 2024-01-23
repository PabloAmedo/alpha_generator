# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 18:11:46 2023

@author: diego
"""

import numpy as np

#<n_{e, cl}> CALCULATION

W_i = 26.4                    #eV
rho_Ar = 1.784e-3             #g/cm3
#next quantities taken from Clusters_calculation.py
dEdx_4GeV = 2.697             #MeV cm2/g
Ncl_4GeV = 55.16              #clusters/cm

n_prom = ( ( dEdx_4GeV ) * rho_Ar ) / ( Ncl_4GeV * W_i * (1e-6) )

print('<n_e,cl>\t=\t', n_prom)
#------------------------------------------------------------------------------

#n_max CALCULATION --> This has to fit the  n_max

n = 1e6
num = 0
denom = 0

for i in range(1,int(n)):
    num = num + 1 / i
    denom = denom + ( 1 / i**2 )
    
    if num/denom > n_prom:
        print('n_max=', i, '\t->\t', num/denom)
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

