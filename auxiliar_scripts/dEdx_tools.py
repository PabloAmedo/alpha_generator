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

def Momentum(energy_list, mass):
    p_list = []
    
    if type(energy_list) == list:
        for energy in energy_list:
            p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        return np.array(p_list)
    else:
        p = np.sqrt((energy_list + mass)**2 - mass**2)
        return p

def load_cd(path):
    data_         =       np.loadtxt( path, skiprows=2, delimiter=' ')
    pm_S_, dNdx_S_ =       np.split(data_, 2, axis = 1)                        #pm  ->  p/m ; 
    pm_S_ = pm_S_
    pm_S = [] ; dNdx_S = []
    for a, b in zip(pm_S_, dNdx_S_):
        pm_S.append(float(a))
        dNdx_S.append(float(b))
        
    return np.array(pm_S), np.array(dNdx_S)



def Get_dNdx(Energy, mass,  Wi = 26.4, Pressure = 10, path = 'alpha_generator/data/cluster_densities/', gas = 'ArCF4', particle = 'muon', pure_argon = False):
    """
    Get the number of ionization clusters per unit length (cl/cm) for given:
        
        - Ekin / momentum
        - Particle (mass)
        - Ionization medium (gas and proportion as Name_XY)
        - Pressure
        
    (!!!) We only have data for ArCF4_99-1. In case we want to get results for 
    other gases (PureAr for instance) we simply escale the number of clusters 
    for a known factor (Santovetti et al. work with calculated dNdx for 
    different noble gases). 
        - For PureAr we can set pure_argon argument = True
    
    """
    try:
        
        fullpath_mu = path + gas + 'cd_' +  'muon' + '.txt'
        fullpath_p = path + gas + 'cd_' +  'proton' + '.txt'
    
        p_mu, dNdx_mu = load_cd(fullpath_mu)                                    #Data for higher momentum
        p_p, dNdx_p = load_cd(fullpath_p)                                       #Data for lower momentum
    except FileNotFoundError:
        
        fullpath_mu = 'data/cluster_densities/' + gas + 'cd_' +  'muon' + '.txt'
        fullpath_p = 'data/cluster_densities/' + gas + 'cd_' +  'proton' + '.txt'
    
        p_mu, dNdx_mu = load_cd(fullpath_mu)                                    #Data for higher momentum
        p_p, dNdx_p = load_cd(fullpath_p)                                       #Data for lower momentum
    
    pm_mu = p_mu / 105.66e6
    pm_p = p_p / 938.272e6
    
    #P/m calculation
    pm = Momentum(Energy, mass) / mass
    #print('pm:\t',pm)
    
    if pm <= 1:
        #print('quad')
        extrapolator = interp1d(pm_p, dNdx_p, kind = 'quadratic', fill_value = 'extrapolate')  
        dNdx_extrapolated = extrapolator(pm)
        if pure_argon:
            dNdx_extrapolated = dNdx_extrapolated / 1.0561177042126448             #scale factor from Ar/CF4 (99/1) to Pure Ar
    else:
        #print('linear')
        data_interpolator = interp1d(pm_mu, dNdx_mu, kind = 'linear', fill_value = 'extrapolate') #slinear
        dNdx_extrapolated = data_interpolator(pm)
        if pure_argon:
            dNdx_extrapolated = dNdx_extrapolated / 1.0561177042126448             #scale factor from Ar/CF4 (99/1) to Pure Ar
    #print('dN:\t',dNdx_extrapolated / Pressure)
    return dNdx_extrapolated / Pressure

"""
def Momentum(energy_list, mass):
    p_list = []
    
    if type(energy_list) == list:
        for energy in energy_list:
            p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        return np.array(p_list)
    else:
        p = np.sqrt((energy_list + mass)**2 - mass**2)
        return p
"""

# PARAMETRIZATIONS OF THE N CLUSTERS PER UNIT LENGTH FOR DIFFERENT GASES BASED ON FISCHLE ET AL. WORK
# Define the 1/n^2 parametrizarions from data
ClusterParametrizationAr = lambda n : 0.216 / n**2 
ClusterParametrizationCH4 = lambda n : 0.119 / n**2 
ClusterParametrizationGeneral = lambda n : 1 / n**2 
