# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:15:26 2024

@author: usuario
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

os.chdir('../') #!!!!!!!!!!

def load_cd(path):
    data_S         =       np.loadtxt( path, skiprows=2, delimiter=' ')
    pm_S_, dNdx_S_ =       np.split(data_S, 2, axis = 1)                       #pm  ->  p/m ; 
    pm_S_ = pm_S_ * 1e-8
    pm_S = [] ; dNdx_S = []
    for a, b in zip(pm_S_, dNdx_S_):
        pm_S.append(float(a))
        dNdx_S.append(float(b))
        
    return pm_S, dNdx_S

def dNdx(Energy, mass,  Wi = 26.4, data = None):
    
    #Load data scanned from paper
    #data_S         =       np.loadtxt('alpha_generator/data/' + data, delimiter=';')  #_S  ->  ref to Santovetti's data
    print(os.getcwd())
    pm_mu, dNdx_mu = load_cd('data/cluster_densities/ArCF4cd_muon.txt')
    pm_pi, dNdx_pi = load_cd('data/cluster_densities/ArCF4cd_pion.txt')
    pm_k, dNdx_k = load_cd('data/cluster_densities/ArCF4cd_kaon.txt')
    pm_p, dNdx_p = load_cd('data/cluster_densities/ArCF4cd_proton.txt')
    
    
    #Momentum calculation
    momentum = np.sqrt((Energy + mass)**2 - mass**2)                           #mass units
    pm = momentum / mass

    muon_interpolator = interp1d(pm_mu, dNdx_mu, kind='linear', fill_value='extrapolate') #slinear
    extrapolator = interp1d(pm_mu, dNdx_mu, kind='zero', fill_value='extrapolate')
    
    data_interpolator = interp1d(pm_mu, dNdx_mu, kind='linear', fill_value='extrapolate') #slinear
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
