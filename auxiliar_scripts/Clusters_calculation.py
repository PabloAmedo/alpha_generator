# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:42:22 2023

@author: diego
"""

import numpy as np
from scipy.interpolate import interp1d

# CALCULATION OF CLUSTERS/MM 
data_path = '../data/'

def dNdx(Energy, mass,  Wi = 26.4):
    
    """
    The following calculation uses the results from Santovetti et al. "Primary 
    ionization and energy loss calculation for helium, neon, argon and krypton" 
    where curves for the number of clusters per cm are given for different 
    gasses.
    """
    #Load data scanned from paper
    data_S       =       np.loadtxt('data/Ncl_Ar_good.csv', delimiter=';')  #_S  ->  ref to Santovetti's data
    pm_S_, dNdx_S_ =       np.split(data_S, 2, axis = 1)                         #pm  ->  p/m ; 
    pm_S = [] ; dNdx_S = []
    for a, b in zip(pm_S_, dNdx_S_):
        pm_S.append(float(a))
        dNdx_S.append(float(b))
        
    #Momentum calculation
    momentum = np.sqrt((Energy + mass)**2 - mass**2)                           #mass units
    pm = momentum / mass
    
    f_extrapolate = interp1d(pm_S, dNdx_S, kind='zero', fill_value='extrapolate')
    dNdx_extrapolated = f_extrapolate(pm)
    
    #TESTING
    #print('p\t=\t', momentum)
    #print('p/m\t=\t', pm)
    #print('N/cm\t=\t', dNdx_extrapolated)
   
    return dNdx_extrapolated

