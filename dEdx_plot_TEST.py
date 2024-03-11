# -*- coding: utf-8 -*-
"""
Created on Fri Feb  28 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from Alpha_track_simulator import *
from scipy.optimize import curve_fit

from scipy.ndimage import gaussian_filter


#INPUTS
n_tracks = 1
low_E = np.linspace(0.01, 0.10, 25)
high_E = np.linspace(0.1, 1000, 25)
energy_list = list(low_E) + list(high_E)  #[0.01, 0.1, 1, 10, 100]   #(MeV)
dimensions = [250,250,250]          #We are considering 2.5 (m) length 
pxs = 100                           #n Pixels
masses = [105.66, 139.57, 493.7, 938.27]                   #(MeV/c2)  [muon, pi, k, proton]



#lcp_list = []                      #lcp stands for: Light Charged Particle (they are likely to be muons, pions, kaons and protons?)
dEdx = []
momentum = []

for energy in energy_list:
    lcp_list = []  
    particle_gen = muon_generator(energy = energy, geometry = dimensions)      #First generate the lcp generator object
    #Particle generation
    particle_gen.produce_muon(n = n_tracks, store = lcp_list, gas= 'Argon')    #generate lcp's obj stored in muons list
    
    #Iteration over each track
    for lcp in lcp_list:
        
        track_length = lcp.track_length
        n_electrons = lcp.n_electrons
        for mass in masses:
            p = np.sqrt((energy + mass)**2 - mass**2)    #MeV/c
            E_loss = n_electrons * 26.4e-6          #MeV
            momentum.append(p)
            dEdx.append(E_loss/track_length)


plt.figure()
plt.plot(momentum, dEdx, 'kx')
plt.xscale('log')
plt.legend()
plt.grid()
     
        