# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:42:22 2023

@author: diego
"""

import numpy as np
from scipy.interpolate import interp1d
import os
os.chdir('../')
# CALCULATION OF CLUSTERS/MM 

def dNdx(Energy, mass,  Wi = 26.4, data = None):
    
    """
    Function to calculate the Cluster Density for a given energy range and mass 
    (particle).
    You can specify which data you want to use for the calculation:
        - Ncl_Ar_good.csv : uses the results from Santovetti et al. 
        "Primary ionization and energy loss calculation for helium, neon, argon 
        and krypton" 
        
        - cluster_densities/[choose]: data from HEED
    """
    #Load data scanned from paper
    #data_S         =       np.loadtxt('alpha_generator/data/' + data, delimiter=';')  #_S  ->  ref to Santovetti's data
    print(os.getcwd())
    data_S         =       np.loadtxt('data/cluster_densities/ArCF4cd_muon.txt', skiprows=2, delimiter=' ')
    pm_S_, dNdx_S_ =       np.split(data_S, 2, axis = 1)                       #pm  ->  p/m ; 
    pm_S_ = pm_S_ * 1e-8
    pm_S = [] ; dNdx_S = []
    for a, b in zip(pm_S_, dNdx_S_):
        pm_S.append(float(a))
        dNdx_S.append(float(b))
        
    #Momentum calculation
    momentum = np.sqrt((Energy + mass)**2 - mass**2)                           #mass units
    pm = momentum / mass
    
    data_interpolator = interp1d(pm_S, dNdx_S, kind='linear', fill_value='extrapolate') #slinear
    extrapolator = interp1d(pm_S, dNdx_S, kind='zero', fill_value='extrapolate')
    pmDataRange = pm[pm < 1000 ]
    pmOverDataRange = pm[pm >= 1000 ]
    
    if pm < 500 and mass != 0.511:
        dNdx_extrapolated = data_interpolator(pm)
    elif pm < 500 and mass == 0.511:
        dNdx_extrapolated = extrapolator(pm)
    else:
        dNdx_extrapolated = extrapolator(pm)
    
    #TESTING
    #print('p\t=\t', momentum)
    #print('p/m\t=\t', pm)
    print('N/cm\t=\t', dNdx_extrapolated/10)
   
    return dNdx_extrapolated/10

def Momentum(energy_list, mass):
    p_list = []
    for energy in energy_list:
        p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)