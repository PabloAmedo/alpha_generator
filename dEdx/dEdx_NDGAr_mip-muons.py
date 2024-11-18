# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:09:49 2024

@author: diego
"""

import os
os.chdir('../')

from Alpha_track_simulator import *
from scipy.optimize import curve_fit
from scipy.stats import trim_mean
from matplotlib.colors import Normalize
import time

st = time.time()

#Avoid warning in console
import warnings
from scipy.optimize import OptimizeWarning
#from matplotlib.colors import DivergingNorm
warnings.filterwarnings("ignore",  category=OptimizeWarning)

#INPUTS========================================================================
gas             =       'ArCF4_99-1'
n_cl_cm         =       False #29.6933                                         #from HEED simulations
       
energy          =       242             #MeV --> mip      
dimensions      =       [500,500,500]   #cm --> NDGAr 
mass            =       105.66          #MeV/c2
Pressure        =       10              #bar
sampling_size   =       0.2             #cm

dE              =       []
ne              =       []                  

n_particles     =       2000
e_cut           =       10000
tracks_list = []
# =============================================================================
#Gas
ArCF4 = Gas(gas = gas, L_drift = 250, Pressure = Pressure)
#Beam
beam = muon_generator(energy = energy, geometry = dimensions, gas= gas, mass = mass, pressure = Pressure)

for i in range(n_particles):
    tracks_list = []
    beam.produce_muon(n = 1, store = tracks_list, line=True, e_cut = e_cut)
    ne.append(tracks_list[0].n_electrons)

ne = np.array(ne)                               #total number of electrons generated per track
dE = ne * ArCF4.Wi  *1e-6                       #MeV
print('<dE> = {:.3f} (MeV)'.format(dE.mean()))

plt.figure()
plt.title('Ecut = {:.2f} (MeV)'.format(e_cut * ArCF4.Wi*1e-6), fontsize = 15)
plt.hist(dE, bins = 150, histtype = 'step')
plt.plot([], [], 'w-', label = 'sim events: {:.2e}'.format(n_particles))
plt.xlabel('Energy Loss (MeV)', fontsize = 15)
plt.legend()

end = time.time()
print('Time elapsed:\t{:.2f} s'.format(end - st))












