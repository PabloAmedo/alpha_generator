# -*- coding: utf-8 -*-
"""
Created on Fri Feb  28 2024

@author: diego
"""
import os
os.chdir('../')

import numpy as np
import matplotlib.pyplot as plt
from Alpha_track_simulator import *
from scipy.optimize import curve_fit
from scipy.stats import trim_mean
from matplotlib.colors import Normalize
import time

#Avoid warning in console
import warnings
from scipy.optimize import OptimizeWarning
#from matplotlib.colors import DivergingNorm
warnings.filterwarnings("ignore",  category=OptimizeWarning)

start = time.time()
print('Running...')

#INPUTS========================================================================
gas             =       'ArCF4-99/1_01mm'
energy_list     =       np.logspace(1, 5, 10000)                               #(MeV)
dimensions      =       [250,250,0]                                            #We are considering 2.5 (m) length 
masses          =       [105.66, 139.57, 493.7, 938.27]                        #(MeV/c2)  [muon, pi, k, proton]
Pressure        =       10                                                     #bar
sampling_size   =       0.2                                                    #cm (2 mm)

dEdx            =       []
momentum        =       []
binsx           =       int(150*2)
binsy           =       int(75*2)                                                   #NOT SURE ABOUT THIS ONE

n_particles     =       250000

def Momentum(energy_list, mass):
    p_list = []
    for energy in energy_list:
        p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)

"""
--> How to make it FASTER? <--

PRODUCE_MUONS RETURNS THE CLUSTERS' POSITION, SO IT MIGHT BE EASIER THAT I'VE 
BEEN DOING.

-I think it's not going to work because it returns a list(array) for ALL the 
tracks simulated at once. (Take a look to this).

-The energy list selection for every particle could be improved (I think this 
is taking too much time).
"""

energy_list_muon = energy_list
energy_list_pi = energy_list[(energy_list > 39)]
energy_list_k = energy_list[(energy_list > 136)]
energy_list_p = energy_list[(energy_list > 260)]

p_mu = Momentum(energy_list_muon, masses[0]) 
p_pi = Momentum(energy_list_pi, masses[1]) 
p_k = Momentum(energy_list_k, masses[2]) 
p_p = Momentum(energy_list_p, masses[3]) 

p_resolution_mu = np.loadtxt('data/p_resolution/p_muon_res', comments='#') / 100
p_resolution_pi = np.loadtxt('data/p_resolution/p_pion_res', comments='#') / 100
p_resolution_k = np.loadtxt('data/p_resolution/p_kaon_res', comments='#') / 100
p_resolution_p = np.loadtxt('data/p_resolution/p_proton_res', comments='#') / 100

momentum_muon = np.array((p_mu, p_resolution_mu, energy_list_muon)).T
momentum_pion = np.array((p_pi, p_resolution_pi, energy_list_pi)).T
momentum_k = np.array((p_k, p_resolution_k, energy_list_k)).T
momentum_p = np.array((p_p, p_resolution_p, energy_list_p)).T
#CALCULATIONS==================================================================

dedx_truncated = []

for i in range(n_particles):
    mass    = np.random.choice(masses)
    
    if mass == masses[0]:    
        indx = np.random.choice(len(energy_list_muon))
        p, sigmaP, energy = momentum_muon[indx]
    elif mass == masses[1]:
        indx = np.random.choice(len(energy_list_pi))
        p, sigmaP, energy = momentum_pion[indx]
    elif mass == masses[2]:
        indx = np.random.choice(len(energy_list_k))
        p, sigmaP, energy = momentum_k[indx]
    else:
        indx = np.random.choice(len(energy_list_p))
        p, sigmaP, energy = momentum_p[indx]
        
    particle_gen = muon_generator(energy = energy, geometry = dimensions, gas= gas, mass = mass, pressure = Pressure)
    tracks_list = []
    particle_gen.produce_muon(n = 1, store = tracks_list, line=True )
    
    #Iteration over each track
    for track in tracks_list:
        track.fill()                                                   #Maybe this is taking too much time (?)
        x_position, _ = np.split(track.electron_positions, 2, axis = 1)#Getting the x_position (along the track) 
        track_length = track.track_length
        n_electrons = track.n_electrons
        i = sampling_size
        electrons_per_sample = []                                      #Initializing the list for every track!!
        while i <= track_length: #This could be done with Pandas??
            n_electrons_sample = len(x_position[(x_position > i) & (x_position < (i + sampling_size))])
            electrons_per_sample.append(n_electrons_sample)
            i = i + sampling_size
        
        electrons_per_sample = np.array(electrons_per_sample)
        
        E_loss = electrons_per_sample * 26.4e-3                                #keV
        momentum.append(np.random.normal(loc = p, scale = sigmaP))
        dEdx.append(np.mean(E_loss) / sampling_size)                           #keV/cm
        dedx_truncated.append(trim_mean(E_loss, 0.2) / sampling_size)
        print('E',track.energy)
        print('p/m', p / mass)

#PLOTS / OUTS==================================================================
momentum = np.array(momentum)
dedx_truncated = np.array(dedx_truncated)

"""
#SCATTER PLOT
plt.figure()
plt.title('Sample {} (cm) // Pressure {} (bar)'.format(sampling_size, Pressure), fontsize = 20)
plt.plot(momentum / 1e3, dedx_truncated, 'kx')
plt.xscale('log')
plt.xlabel('p (GeV/c)',     fontsize = 13)
plt.ylabel('dEdx (keV/cm)', fontsize = 13)
plt.grid()
"""
#Define binning for dedx v p plot
x_space = np.logspace(0.1,6, binsx)                                            #FIXME -- this is on MeV/c
y_space = np.linspace(min(dedx_truncated) - 0.1 * min(dedx_truncated), max(dedx_truncated) + 0.1 * max(dedx_truncated), binsy)
#norm = DivergingNorm(vmin=dedx_truncated.min(), vcenter=0, vmax=dedx_truncated.max())

plt.figure()
plt.hist2d(momentum, dedx_truncated, bins=(x_space, y_space), cmin=1,cmap = "jet" )
plt.xscale('log')
plt.colorbar()

#MAYBE I SHOULD NORMALIZE THE VALUES FOR THE CMAP THINGY --> You can edit it in the plot gui

end = time.time()
execution_time = end-start
print('Time elapsed:\t', execution_time)

