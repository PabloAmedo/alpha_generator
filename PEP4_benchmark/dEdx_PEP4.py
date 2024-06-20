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


def Momentum(energy_list, mass):
    p_list = []
    for energy in energy_list:
        p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)
os.chdir('alpha_generator/')

#INPUTS========================================================================
gas             =       'PEP4'
n_cl_cm         =       False #29.6933                                         #from HEED simulations
name            =       'PEP4_clust'
file_tosave     =       open('PEP4_benchmark/simulated_data/'+ name +'.txt', 'x');         
file_tosave.write('\nINPUTS\n'+'='*50+'\n')  
energy_list     =       np.logspace(1, 5, 10000) ;          file_tosave.write('Energy range:\t{} - {} ({}) (MeV)\n'.format(energy_list[0], energy_list[-1], len(energy_list)))
dimensions      =       [100,100,0];                        file_tosave.write('Dimensions:\t{} (cm)\n'.format(dimensions)) #We are considering 2.5 (m) length 
masses          =       [0.511, 105.66, 139.57, 493.7, 938.27];    file_tosave.write('Masses:\t\t{} (MeV)\n'.format(masses))  #(MeV/c2)  [muon, pi, k, proton]
Pressure        =       8.5;                                 file_tosave.write('Pressure:\t{} (bar)\n'.format(Pressure))#bar
sampling_size   =       0.4;                                file_tosave.write('Sampling size:\t{} (cm)\n\n'.format(sampling_size)) #cm (2 mm) ; #4 mm in PEP-4

dEdx            =       []
momentum        =       []
binsx           =       int(150*2)
binsy           =       int(75*2);                          file_tosave.write('Bins(x,y):\t{}\n'.format([binsx, binsy]))

n_particles     =       30000;                                file_tosave.write('SIMULATED EVENTS:\t{}\n'.format(n_particles))




file_tosave.write('\n' + '='*50 + '\n')


#Cuts in energy to match the first p/m for every particle --> avoid interpolation at low p
energy_list_muon = energy_list
energy_list_pi = energy_list[(energy_list > 39)]
energy_list_k = energy_list[(energy_list > 136)]
energy_list_p = energy_list[(energy_list > 260)]

p_e = Momentum(energy_list, masses[0])
p_mu = Momentum(energy_list_muon, masses[1]) 
p_pi = Momentum(energy_list_pi, masses[2]) 
p_k = Momentum(energy_list_k, masses[3]) 
p_p = Momentum(energy_list_p, masses[4]) 

p_resolution_e = np.loadtxt('data/p_resolution/p_electron_res', comments='#') / 100
p_resolution_mu = np.loadtxt('data/p_resolution/p_muon_res', comments='#') / 100
p_resolution_pi = np.loadtxt('data/p_resolution/p_pion_res', comments='#') / 100
p_resolution_k = np.loadtxt('data/p_resolution/p_kaon_res', comments='#') / 100
p_resolution_p = np.loadtxt('data/p_resolution/p_proton_res', comments='#') / 100

momentum_e = np.array((p_e, p_resolution_e, energy_list)).T
momentum_muon = np.array((p_mu, p_resolution_mu, energy_list_muon)).T
momentum_pion = np.array((p_pi, p_resolution_pi, energy_list_pi)).T
momentum_k = np.array((p_k, p_resolution_k, energy_list_k)).T
momentum_p = np.array((p_p, p_resolution_p, energy_list_p)).T
#CALCULATIONS==================================================================

dedx_truncated = []     #truncated mean
file_tosave.write('Mass\tP_true\tP\tdEdx\tdEdx_trunc\n')

for i in range(n_particles):
    mass    = np.random.choice(masses)
    
    if mass == masses[1]:    
        indx = np.random.choice(len(energy_list_muon))
        p_true, sigmaP, energy = momentum_muon[indx]
    elif mass == masses[2]:
        indx = np.random.choice(len(energy_list_pi))
        p_true, sigmaP, energy = momentum_pion[indx]
    elif mass == masses[3]:
        indx = np.random.choice(len(energy_list_k))
        p_true, sigmaP, energy = momentum_k[indx]
    elif mass == masses[4]:
        indx = np.random.choice(len(energy_list_p))
        p_true, sigmaP, energy = momentum_p[indx]
    else:
        indx = np.random.choice(len(energy_list))
        p_true, sigmaP, energy = momentum_e[indx]
    
    print('mass: ',mass)
    particle_gen = muon_generator(energy = energy, geometry = dimensions, gas= gas, mass = mass, pressure = Pressure)
    tracks_list = []
    particle_gen.produce_muon(n = 1, store = tracks_list, line=True, n_cl_cm_in = n_cl_cm )
    
    #Iteration over each track
    for track in tracks_list:
        track.fill()                                                           #Maybe this is taking too much time (?)
        x_position, _ = np.split(track.electron_positions, 2, axis = 1)        #Getting the x_position (along the track) 
        track_length = track.track_length
        n_electrons = track.n_electrons
        i = sampling_size
        electrons_per_sample = []                                      #Initializing the list for every track!!
        while i <= track_length: #This could be done with Pandas??
            n_electrons_sample = len(x_position[(x_position > i) & (x_position < (i + sampling_size))])
            electrons_per_sample.append(n_electrons_sample)
            i = i + sampling_size
        
        electrons_per_sample = np.array(electrons_per_sample)
        
        E_loss = electrons_per_sample * 26.4e-3   
        p = np.random.normal(loc = p_true, scale = np.sqrt((sigmaP * p_true)**2 + 0.01*(p_true)**2)) #momentum res + 1% de multiple scattering
        momentum.append(p)
        dEdx_ = np.mean(E_loss) / sampling_size
        dEdx.append(dEdx_)                                                     #keV/cm
        dEdx_trunc = trim_mean(E_loss, 0.2) / sampling_size                    #80% !
        dedx_truncated.append(dEdx_trunc)
        
        file_tosave.write('{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(mass, p_true, p, dEdx_, dEdx_trunc))

#PLOTS / OUTS==================================================================
momentum = np.array(momentum)
dedx_truncated = np.array(dedx_truncated)

#Define binning for dedx v p plot
x_space = np.logspace(0.1,6, binsx)                                            #FIXME -- this is on MeV/c
y_space = np.linspace(min(dedx_truncated) - 0.1 * min(dedx_truncated), max(dedx_truncated) + 0.1 * max(dedx_truncated), binsy)

plt.figure()
plt.hist2d(momentum, dedx_truncated, bins=(x_space, y_space), cmin=1,cmap = "jet" )
plt.xscale('log')
plt.colorbar()

#MAYBE I SHOULD NORMALIZE THE VALUES FOR THE CMAP THINGY --> You can edit it in the plot gui

end = time.time()
execution_time = end-start
print('Time elapsed:\t', execution_time)
file_tosave.close()
