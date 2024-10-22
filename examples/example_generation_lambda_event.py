# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:17:39 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('../')
from Alpha_track_simulator import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import time
import matplotlib

t0 = time.time()

#INPUTS========================================================================
p3D         =   1
P           =   10                                                              #bar                                                      #E = 4000 [MeV] (CR)
dimensions  =   [250,250,500]

mass_mu     =   105.66;     E_mu =  2300        #MeV/c2; MeV
mass_pi     =   134.98;     E_pi =  103.18      #MeV/c2; MeV
mass_p      =   938.27;     E_p  =  115.70      #MeV/c2; MeV

gas         =   'ArCF4_99-1'
e_cut       =   10000

px_size     =   0.2

sigma_diff  =   0.21
sigma_PSF   =   0

bins        =   [int(dimensions[0] / px_size), int(dimensions[1] / px_size)]
x_range = [0, (px_size * bins[0])]
y_range = [0, (px_size * bins[1])]

muons       =   []
pions       =   []
protons     =   []

#==============================================================================
#==============================================================================

argon = Gas(gas = 'ArCF4', L_drift = abs(dimensions[1]/2), Pressure = P)

muon_gen = muon_generator(energy = E_mu, geometry = dimensions, gas = gas, pressure = P)
pion_gen = muon_generator(energy = E_pi, geometry = dimensions, gas = gas, pressure = P)
proton_gen = muon_generator(energy = E_p, geometry = dimensions, gas = gas, pressure = P)

muon_gen.produce_muon(n = 1, store = muons, e_cut = e_cut,
                  position_in = [25, int(np.random.rand() * dimensions[1] / 2), int(np.random.rand() * dimensions[2] / 2)])
pion_gen.produce_muon(n = 1, store = pions, e_cut = e_cut,
                  position_in = [70, 121, 150], ath_in = np.pi/2, phi_in = -0.36)#, length = 15)
proton_gen.produce_muon(n = 1, store = protons, e_cut = e_cut,
                  position_in = [70, 70, 150], ath_in = np.pi/2, phi_in = 0.146)#, length = 15)

#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff = sigma_diff, sigma_PSF=0)
diff_handler.diffuse(muons)
diff_handler.diffuse(pions)
diff_handler.diffuse(protons)

#Filling tracks with e-
for muon in muons: muon.fill(diff = True)
for pion in pions: pion.fill(diff = True)
for proton in protons: proton.fill(diff = True)

#Noise object
noise=Noise(10)
#Plot the final tracks
image2d_mu = Image_2D(track_list = muons, hist_args={"bins":bins, "range": [x_range, y_range]}, pixel_size = px_size)
image2d_pi = Image_2D(track_list = pions, hist_args={"bins":bins, "range": [x_range, y_range]}, pixel_size = px_size)
image2d_p = Image_2D(track_list = protons, hist_args={"bins":bins, "range": [x_range, y_range]}, pixel_size = px_size)
plt.close('all')

#Getting histogram and edges to plot with white background
hist_mu = image2d_mu.Hist2D_e
hist_x_mu = image2d_mu.x_edges
hist_y_mu = image2d_mu.y_edges

hist_pi = image2d_pi.Hist2D_e
hist_p = image2d_p.Hist2D_e

hist = hist_mu + hist_pi + hist_p
hist[hist == 0] = 'NaN'                                                        #pcolormesh plots nans as white

plt.figure()
plt.pcolormesh(hist_x_mu, hist_y_mu, hist, cmap = 'jet')
plt.colorbar()
plt.xlim([0, dimensions[0]])
plt.ylim([0, dimensions[1]])
plt.xlabel('X (cm)', fontsize = 15)
plt.ylabel('Y (cm)', fontsize = 15)

print('px size:\t({:.3f} x {:.3f}) (cm x cm)'.format(np.mean(np.diff(hist_x_mu)), np.mean(np.diff(hist_y_mu))))
#3D-PLOT ======================================================================
if p3D == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = []
    y = []
    z = []
    colors = colors = list(matplotlib.colors.XKCD_COLORS)
    for i in range(len(muons)):
        x_mu = muons[i].electron_positions_diff[:,0]
        y_mu = muons[i].electron_positions_diff[:,1]
        z_mu = muons[i].electron_positions_diff[:,2]
        
        x_pi = pions[i].electron_positions_diff[:,0]
        y_pi = pions[i].electron_positions_diff[:,1]
        z_pi = pions[i].electron_positions_diff[:,2]
        
        x_p = protons[i].electron_positions_diff[:,0]
        y_p = protons[i].electron_positions_diff[:,1]
        z_p = protons[i].electron_positions_diff[:,2]
    
        ax.plot(x_mu, y_mu, z_mu,'r.')
        ax.plot(x_pi, y_pi, z_pi,'g.')
        ax.plot(x_p, y_p, z_p,'b.')
    cube_size = dimensions[0]
    
    for i in range(2):
        for j in range(2):
            ax.plot([0, dimensions[0]], [i*dimensions[1], i*dimensions[1]], [j*dimensions[2], j*dimensions[2]], color='k')
            ax.plot([i*dimensions[0], i*dimensions[0]], [0, dimensions[1]], [j*dimensions[2], j*dimensions[2]], color='k')
            ax.plot([i*dimensions[0], i*dimensions[0]], [j*dimensions[1], j*dimensions[1]], [0, dimensions[2]], color='k')
            
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')
    max_dim = max(dimensions)
    ax.set_xlim([0 - 3 , max_dim + 3])
    ax.set_ylim([0 - 3 , max_dim + 3])
    ax.set_zlim([0 - 3, max_dim + 3])
    ax.grid(False)
    plt.show()

tf = time.time()

print('Time elapsed:\t{} s'.format(tf - t0))











