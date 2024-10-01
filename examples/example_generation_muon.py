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

def Momentum(energy, mass):
    p_list = []
    p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)

#INPUTS========================================================================
n_tracks    =   1
P           =   1.5                                                          #bar
#y0          =   10
#theta0      =   50
E           =   2500                                                          #E = 4000 [MeV] (CR)
dimensions  =   [14,14,14]
mass        =   105.66
gas         =   'ArCF4_99-1_01mm'
e_cut       =   10000

px          =   13e-3 
rebin       =   4 
mag         =   11 
px_size     =   (px * rebin * mag)*1e-1

sigma_diff  =   0.17    #!!!!!!!!!!!!!!!!!! NOT SURE
sigma_PSF   =   0
line        =   False

bins        =   [int(dimensions[0] / px_size), int(dimensions[1] / px_size)]
x_range = [0, (px_size * bins[0])]
y_range = [0, (px_size * bins[1])]

muons       =   [] 

#==============================================================================
#==============================================================================

#Muons generation
muon = muon_generator(energy = E, geometry = dimensions, gas = gas, pressure = P)            #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons, e_cut = e_cut, line = line)     #generate muon's obj stored in muons list

#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff= sigma_diff ,sigma_PSF=sigma_PSF)
diff_handler.diffuse(muons)
#Filling tracks with e-
for muon in muons: muon.fill()

#Noise object
noise=Noise(10)
#Plot the final tracks
image2d=Image_2D(track_list = muons, hist_args={"bins":bins, "range": [x_range, y_range]}, pixel_size = px_size)
image2d.track_plot()
image2d.plot_hist(noise_object=noise)
image2d.plot_x()

plt.close('all')

#Getting histogram and edges to plot with white background
hist = image2d.Hist2D_e
hist_x = image2d.x_edges
hist_y = image2d.y_edges

hist[hist == 0] = 'NaN'                                                        #pcolormesh plots nans as white

plt.figure()
plt.pcolormesh(hist_x, hist_y, hist, cmap = 'jet')
plt.colorbar()
plt.xlim([0, dimensions[0]])
plt.ylim([0, dimensions[1]])
plt.xlabel('X (cm)', fontsize = 15)
plt.ylabel('Y (cm)', fontsize = 15)

print('px size:\t({:.3f} x {:.3f}) (cm x cm)'.format(np.mean(np.diff(hist_x)), np.mean(np.diff(hist_y))))

plt.figure()
plt.imshow(image2d.Hist2D)