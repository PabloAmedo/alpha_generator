# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:16:43 2023

@author: diego
"""
import os
os.chdir('../')
from Alpha_track_simulator import *
import scipy.special as spc


plt.close('all')
print('Running...\n')

#INPUTS========================================================================
n_tracks    =   1
P           =   10                                                             #bar
line        =   False
E           =   4000                                                           #E = 4000 [MeV] (CR)
dimensions  =   [50,50,50]
mass        =   105.66
gas         =   'Argon'
e_cut       =   1000

bins        =   100
sigma_diff  =   0.18
sigma_PSF   =   0

muons       =   [] 

#==============================================================================
#==============================================================================

#Muons generation
muon = muon_generator(energy = E, geometry = dimensions, gas = gas, pressure = P)            #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons, line = line, e_cut = e_cut)     #generate muon's obj stored in muons list

#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff= sigma_diff ,sigma_PSF=sigma_PSF)
diff_handler.diffuse(muons)
#Filling tracks with e-
for muon in muons: muon.fill()

#Noise object
noise=Noise(10)
#Plot the final tracks
image2d=Image_2D(track_list = muons, hist_args={"bins":100})
image2d.track_plot()
image2d.plot_hist(noise_object=noise)
image2d.plot_x()
