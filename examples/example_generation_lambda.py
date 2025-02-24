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
warnings.filterwarnings("ignore", category = DeprecationWarning)
import time
import matplotlib

t0 = time.time()

#INPUTS========================================================================
p3D         =   1
P           =   1                               #bar
dimensions  =   [30,30,30]

mass_mu     =   105.66;     E_mu =  2300        #MeV/c2; MeV
mass_pi     =   134.98;     E_pi =  103.18      #MeV/c2; MeV
mass_p      =   938.27;     E_p  =  115.70      #MeV/c2; MeV

gas         =   'ArCF4_99-1'
e_cut       =   10000

px_size     =   0.2
bins3d      =   50

sigma_diff  =   0.21
sigma_PSF   =   0

bins        =   [int(dimensions[0] / px_size), int(dimensions[1] / px_size)]
x_range     =   [0, (px_size * bins[0])]
y_range     =   [0, (px_size * bins[1])]

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
                  position_in = [10*0.3, 15*0.3,27*0.3],
                  phi_in = -0.08489, ath_in = np.pi/2)
pion_gen.produce_muon(n = 1, store = pions, e_cut = e_cut,
                  position_in = [35*0.3, 20*0.3, 30*0.3], ath_in =  0.75 *np.pi/2, phi_in = 0.543337) #, length = 15)
proton_gen.produce_muon(n = 1, store = protons, e_cut = e_cut,
                  position_in = [35*0.3, -40*0.3, 30*0.3], ath_in = 0.75 * np.pi/2, phi_in = 1.011203)#, length = 15)

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
image2d_mu = Image_2D(track_list = muons, 
                      hist_args = {"bins":bins, "range": [x_range, y_range]})#, pixel_size = px_size)
image2d_pi = Image_2D(track_list = pions, 
                      hist_args = {"bins":bins, "range": [x_range, y_range]})#, pixel_size = px_size)
image2d_p = Image_2D(track_list = protons, 
                     hist_args = {"bins":bins, "range": [x_range, y_range]})#, pixel_size = px_size)
plt.close('all')

#Getting histogram and edges to plot with white background
hist_x_mu = image2d_mu.x_edges
hist_y_mu = image2d_mu.y_edges

hist_mu = image2d_mu.Hist2D_e
hist_pi = image2d_pi.Hist2D_e
hist_p = image2d_p.Hist2D_e
hist = hist_mu + hist_pi + hist_p
hist[hist == 0] = 'NaN'    

print('px size:\t({:.3f} x {:.3f}) (cm x cm)'.format(np.mean(np.diff(hist_x_mu)), np.mean(np.diff(hist_y_mu))))

# =============================================================================
# =============================================================================
# PLOT
# =============================================================================

# Crear la figura para ambos subplots
fig = plt.figure(figsize=(10, 5))

# First subplot: XY projection as seen by the camera
ax1 = fig.add_subplot(1, 2, 1)
pcolormesh_plot = ax1.pcolormesh(hist_x_mu, hist_y_mu, hist, cmap='hot')
plt.colorbar(pcolormesh_plot, ax=ax1)
ax1.set_xlim([0, dimensions[0]])
ax1.set_ylim([0, dimensions[1]])
ax1.set_xlabel('X (cm)', fontsize=15)
ax1.set_ylabel('Y (cm)', fontsize=15)

#------------------------------------------------------------------------------
#Second subplot: 3D plot
ax2 = fig.add_subplot(1, 2, 2, projection='3d')

all_x = []
all_y = []
all_z = []

#Getting coordinates of all particles in the same array
for i in range(len(muons)):
    #muon
    x_mu = muons[i].electron_positions_diff[:, 0]
    y_mu = muons[i].electron_positions_diff[:, 1]
    z_mu = muons[i].electron_positions_diff[:, 2]
    all_x.append(x_mu)
    all_y.append(y_mu)
    all_z.append(z_mu)
    #pion
    x_pi = pions[i].electron_positions_diff[:, 0]
    y_pi = pions[i].electron_positions_diff[:, 1]
    z_pi = pions[i].electron_positions_diff[:, 2]
    all_x.append(x_pi)
    all_y.append(y_pi)
    all_z.append(z_pi)
    #proton
    x_p = protons[i].electron_positions_diff[:, 0]
    y_p = protons[i].electron_positions_diff[:, 1]
    z_p = protons[i].electron_positions_diff[:, 2]
    all_x.append(x_p)
    all_y.append(y_p)
    all_z.append(z_p)

x_all = np.concatenate(all_x)
y_all = np.concatenate(all_y)
z_all = np.concatenate(all_z)

#3D histogram to get density --> colormap
hist, edges = np.histogramdd(np.array([x_all, y_all, z_all]).T, bins=(bins3d, bins3d, bins3d))  # 50x50x50 binning
density = hist.flatten()

#3D mesh -> squares of hist-bins size
x_bin_centers = 0.5 * (edges[0][:-1] + edges[0][1:])
y_bin_centers = 0.5 * (edges[1][:-1] + edges[1][1:])
z_bin_centers = 0.5 * (edges[2][:-1] + edges[2][1:])
x_mesh, y_mesh, z_mesh = np.meshgrid(x_bin_centers, y_bin_centers, z_bin_centers, indexing='ij')

density_map = np.zeros_like(x_all)

#Getting the ionization density
for i in range(len(x_all)):
    
    bin_x = np.digitize(x_all[i], x_bin_centers) - 1
    bin_y = np.digitize(y_all[i], y_bin_centers) - 1
    bin_z = np.digitize(z_all[i], z_bin_centers) - 1
    
    # Check that the points are inside
    if 0 <= bin_x < len(x_bin_centers) and 0 <= bin_y < len(y_bin_centers) and 0 <= bin_z < len(z_bin_centers):
        density_map[i] = density[bin_x * len(y_bin_centers) * len(z_bin_centers) + 
                                  bin_y * len(z_bin_centers) + bin_z]
        
#Normalization --> we could skip this step
norm_density = (density_map - np.min(density_map)) / (np.max(density_map) - np.min(density_map))

#Plot tracks using the density as cmap-guide
for i in range(len(muons)):
    #muons
    ax2.scatter(muons[i].electron_positions_diff[:, 0],
                muons[i].electron_positions_diff[:, 1],
                muons[i].electron_positions_diff[:, 2],
                c=norm_density[:len(muons[i].electron_positions_diff)], marker='o', cmap='hot', s=1, label='$\\mu^+$' if i == 0 else "")
    
    #pions
    ax2.scatter(pions[i].electron_positions_diff[:, 0],
                pions[i].electron_positions_diff[:, 1],
                pions[i].electron_positions_diff[:, 2],
                c=norm_density[len(muons[i].electron_positions_diff):len(muons[i].electron_positions_diff) + len(pions[i].electron_positions_diff)], marker='o', cmap='hot', s=1, label='$\\pi^-$' if i == 0 else "")
    
    #protons
    ax2.scatter(protons[i].electron_positions_diff[:, 0],
                protons[i].electron_positions_diff[:, 1],
                protons[i].electron_positions_diff[:, 2],
                c=norm_density[len(muons[i].electron_positions_diff) + len(pions[i].electron_positions_diff):], marker='o', cmap='hot', s=1, label='p' if i == 0 else "")

#Plot volume
for i in range(2):
    for j in range(2):
        ax2.plot([0, dimensions[0]], [i * dimensions[1], i * dimensions[1]], [j * dimensions[2], j * dimensions[2]], color='k')
        ax2.plot([i * dimensions[0], i * dimensions[0]], [0, dimensions[1]], [j * dimensions[2], j * dimensions[2]], color='k')
        ax2.plot([i * dimensions[0], i * dimensions[0]], [j * dimensions[1], j * dimensions[1]], [0, dimensions[2]], color='k')

#Labels + limits
ax2.set_xlabel('X (cm)')
ax2.set_ylabel('Y (cm)')
ax2.set_zlabel('Z (cm)')
max_dim = max(dimensions)
ax2.set_xlim([0 - 3, max_dim + 3])
ax2.set_ylim([0 - 3, max_dim + 3])
ax2.set_zlim([0 - 3, max_dim + 3])
ax2.grid(False)

plt.tight_layout()
plt.show()

tf = time.time()

print('Time elapsed:\t{} s'.format(tf - t0))


