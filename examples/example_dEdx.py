# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 16:21:08 2024

@author: diego
"""

import os
os.chdir('../')
from Alpha_track_simulator import *

# INPUTS ======================================================================
n_tracks = 10000
dimensions = [1.5, 1.5, 1.5] #cm
gas = 'ArCH4_93-7'
P = 1 #bar
mass = 139.57 #MeV/c2
energy = 3500
e_cut = 1000000
# =============================================================================
print(os.getcwd())
os.chdir('alpha_generator/')
beam =  muon_generator(energy = energy, geometry = dimensions, gas= gas, mass = mass, pressure = P)
ne = []
for i in range(n_tracks):
    track_list = []
    beam.produce_muon(n = 1, store = track_list, line=True, e_cut = e_cut)
    ne.append(track_list[0].n_electrons)

ne = np.array(ne)
dE = ne * 26.7e-3 #keV

data_Harris = np.loadtxt('data/Harris72.csv', delimiter = ';')

hval, hbin = np.histogram(dE, bins = 50)

plt.figure()
plt.title('Ecut = {:.2f} (MeV)'.format(e_cut * 26.7e-6))
plt.hist(dE, bins = 50, histtype = 'step', label = 'Simulated')
plt.plot(data_Harris[:,0], (data_Harris[:,1]/max(data_Harris[:,1])) * max(hval), label = 'Harris')
plt.legend()
    
    

