# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 16:31:23 2024

@author: diego
"""


import numpy as np
import matplotlib.pyplot as plt
from Alpha_track_simulator import *
from scipy.optimize import curve_fit


plt.close('all')

#INPUTS========================================================================
n_tracks = 1
Emuon = 4000                    #(MeV) // assuming CR as source
muons = []                      #list to store muon objects
dimension = ([5,5,5])           #[x,y,z] (cm) of a rectangular-wise chamber
sigma_diff = 0.015
sigma_PSF = 0
pxs = 100

###############################################################################
#Pixel size calculation
px_size = dimension[0] / pxs       #cm / px


###############################################################################
"""
First of all we have to generate some muons so we use the simulator.
"""
print('Running...')
#Muons generation
muon = muon_generator(energy = Emuon, geometry = dimension)                    #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons, line= True)                                 #generate muon's obj stored in muons list
#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff= sigma_diff ,sigma_PSF=sigma_PSF)
diff_handler.diffuse(muons)
#Filling tracks with e-
for muon in muons: muon.fill()
#Noise object
noise=Noise()
 

bins_x = int(muons[0].x * (1 / px_size))
bins_y = int(abs(muons[0].y - muons[0].y0) * (1 / px_size))

#Plot the final tracks
image2d=Image_2D(track_list = muons, hist_args={"bins":[bins_x, bins_y]})
image2d.track_plot()
image2d.plot_hist(noise_object=noise)
image2d.plot_x()
###############################################################################
e_2d = image2d.Hist2D       #e- positions in the 2d array (detector's pixels (?) )

total_e_col = []            #total electrons in a column of pixels
probab_centroid = []
total_e_col.append(np.sum(e_2d, axis= 1))




