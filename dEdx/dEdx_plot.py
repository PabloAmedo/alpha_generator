# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 17:28:16 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from Alpha_track_simulator import *
from scipy.optimize import curve_fit

#Avoid warning in console
import warnings
from scipy.optimize import OptimizeWarning
warnings.filterwarnings("ignore", category = OptimizeWarning)                  #this avoids the covariance calculation warning 

print('Running...')
"""
NOT WORKING BECAUSE I HAVE TO UPDATE THE NUMBER OF CLUSTERS PER UNIT LENGTH IN 
THE Alpha_track_simulator.py FILE. 
"""
#INPUTS
n_tracks = 1
low_E = np.linspace(0.0001, 100, 5000)
high_E =[]# np.linspace(1, 1000, 500)
energy_list = list(low_E) + list(high_E)  #[0.01, 0.1, 1, 10, 100]   #(MeV)
dimensions = [250,250,250]          #We are considering 2.5 (m) length 
pxs = 100                           #n Pixels
lcp_mass_muon = 105.66                   #(MeV/c2)  MUON 
lcp_mass_pion = 139.57
lcp_mass_k = 493.7


#lcp_list = []                      #lcp stands for: Light Charged Particle (they are likely to be muons, pions, kaons and protons?)
dEdx_muon = []
dEdx_pion = []
dEdx_k = []
momentum_muon = []
momentum_pion = []
momentum_k = []

for energy in energy_list:
    lcp_list = []  
    particle_gen = muon_generator(energy = energy, geometry = dimensions)      #First generate the lcp generator object
    #Particle generation
    particle_gen.produce_muon(n = n_tracks, store = lcp_list, gas= 'ArCH4-90/10', line= True)    #generate lcp's obj stored in muons list
    
    #Iteration over each track
    for lcp in lcp_list:
        
        track_length = lcp.track_length
        n_electrons = lcp.n_electrons
        #MUONS
        p_muon = np.sqrt((energy + lcp_mass_muon)**2 - lcp_mass_muon**2)    #MeV/c
        E_loss_muon = n_electrons * 26.4e-6          #MeV
        momentum_muon.append(p_muon)
        dEdx_muon.append(E_loss_muon/track_length)
        #PIONS
        p_pion = np.sqrt((energy + lcp_mass_pion)**2 - lcp_mass_pion**2)    #MeV/c
        E_loss_pion = n_electrons * 26.4e-6          #MeV
        momentum_pion.append(p_pion)
        dEdx_pion.append(E_loss_pion/track_length)
        #KS
        p_k = np.sqrt((energy + lcp_mass_k)**2 - lcp_mass_k**2)    #MeV/c
        E_loss_k = n_electrons * 26.4e-6          #MeV
        momentum_k.append(p_k)
        dEdx_k.append(E_loss_k/track_length)
    
plt.figure()
plt.plot(momentum_muon, dEdx_muon, 'kx', label = '$\\mu$')
plt.plot(momentum_pion, dEdx_pion, 'bo', label = '$\\pi$')
plt.plot(momentum_k, dEdx_k, 'gv', label ='k')
plt.xscale('log')
plt.legend()
plt.grid()
     
        