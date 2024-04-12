# -*- coding: utf-8 -*-
"""
Created on Fri Feb  28 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from Alpha_track_simulator import *
from scipy.optimize import curve_fit

#Avoid warning in console
import warnings
from scipy.optimize import OptimizeWarning
warnings.filterwarnings("ignore",  category=OptimizeWarning)

print('Running...')

"""
THIS SCRIPT IS SIMULATING TRACKS AS STRAIGHT LINES!!!
"""
#INPUTS========================================================================
n_tracks        =       1
energy_list     =       np.logspace(1,5, 20)                                  #(MeV)
dimensions      =       [250,250,250]                                          #We are considering 2.5 (m) length 
pxs             =       100                                                    #n Pixels
masses          =       [105.66, 139.57, 493.7, 938.27]                        #(MeV/c2)  [muon, pi, k, proton]
bins            =       100 
dEdx            =       []
momentum        =       []
Pressure        =       10                                                     #bar
sampling_sizes  =       [0.2]                                     #cm


"""
TO FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PRODUCE_MUONS RETURNS THE CLUSTERS' POSITION, SO IT MIGHT BE EASIER THAT I'VE 
BEEN DOING. 
"""
#CALCULATIONS==================================================================
for sampling_size in sampling_sizes:
    momentum=[]
    dEdx=[]
    for energy in energy_list:
        for mass in masses:
            particle_gen = muon_generator(energy = energy, geometry = dimensions, gas= 'ArCH4-90/10', mass = mass, pressure = Pressure)  #First generate the lcp generator object
            #Particle generation
            tracks_list = []
            particle_gen.produce_muon(n = n_tracks, store = tracks_list, line=True )    #generate lcp's obj stored in muons list
        
            #Iteration over each track
            for track in tracks_list:
                track.fill()                                                       #Maybe this is taking too much time (?)
                x_position, _ = np.split(track.electron_positions, 2, axis = 1)    #Getting the x_position (along the track) 
                track_length = track.track_length
                n_electrons = track.n_electrons
                i = sampling_size
                electrons_per_sample = []                                          #Initializing the list for every track!!
                while i <= track_length:
                    n_electrons_sample = len(x_position[(x_position > i) & (x_position < (i + sampling_size))])
                    electrons_per_sample.append(n_electrons_sample)
                    i = i + sampling_size
                
                electrons_per_sample = np.array(electrons_per_sample)
                p = np.sqrt((energy + mass)**2 - mass**2)                          #MeV/c
                E_loss = electrons_per_sample * 26.4e-3                            #keV
                momentum.append(p)
                dEdx.append(np.mean(E_loss) / track_length)                        #keV/cm
                print('E',track.energy)
    
    #PLOTS / OUTS==================================================================
    momentum = np.array(momentum)
    
    plt.figure()
    plt.title('Sample {} cm'.format(sampling_size))
    plt.plot(momentum / 1e3, dEdx, 'kx')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    
    
    
    plt.figure()
    plt.title('Sample {} cm'.format(sampling_size))
    hval, hbins, _ = plt.hist(E_loss, bins = bins, histtype='step')
     
        