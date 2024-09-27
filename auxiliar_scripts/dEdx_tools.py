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
    data_         =       np.loadtxt( path, skiprows=2, delimiter=' ')
    pm_S_, dNdx_S_ =       np.split(data_, 2, axis = 1)                        #pm  ->  p/m ; 
    pm_S_ = pm_S_
    pm_S = [] ; dNdx_S = []
    for a, b in zip(pm_S_, dNdx_S_):
        pm_S.append(float(a))
        dNdx_S.append(float(b))
        
    return np.array(pm_S), np.array(dNdx_S)

def dNdx(Energy, mass,  Wi = 26.4, Pressure = 10, path = 'alpha_generator/data/cluster_densities/', file = 'ArCF4cd_', particle = 'muon', ext = '.txt', pure_argon = False):
    print(os.getcwd())
    fullpath = path + file + particle + ext
    
    p_mu, dNdx_mu = load_cd(fullpath)                              #Data for higher p
    p_p, dNdx_p = load_cd(path + file + 'proton' + ext)            #Data for lower p
    
    pm_mu = p_mu / 105.66e6
    pm_p = p_p / 938.272e6
    #Momentum calculation
    momentum = np.sqrt((Energy + mass)**2 - mass**2)                           #mass units
    print('p=', momentum)
    pm = momentum / mass
    print('p/m=', pm)
    
    if pm <= 1:
        extrapolator = interp1d(pm_p, dNdx_p, kind = 'quadratic', fill_value = 'extrapolate')  
        dNdx_extrapolated = extrapolator(pm)
        if pure_argon == True:
            dNdx_extrapolated = dNdx_extrapolated / 1.0561177042126448             #scale factor from Ar/CF4 (99/1) to Pure Ar
    else:
        data_interpolator = interp1d(pm_mu, dNdx_mu, kind = 'linear', fill_value = 'extrapolate') #slinear
        dNdx_extrapolated = data_interpolator(pm)
        if pure_argon == True:
            dNdx_extrapolated = dNdx_extrapolated / 1.0561177042126448             #scale factor from Ar/CF4 (99/1) to Pure Ar

    return dNdx_extrapolated / Pressure

def Momentum(energy_list, mass):
    p_list = []
    
    if type(energy_list) == list:
        for energy in energy_list:
            p_list.append(np.sqrt((energy + mass)**2 - mass**2))
            return np.array(p_list)
    else:
        p = np.sqrt((energy_list + mass)**2 - mass**2)
        return p
