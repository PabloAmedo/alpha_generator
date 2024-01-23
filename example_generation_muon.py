# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:16:43 2023

@author: diego
"""

from Alpha_track_simulator import *
import numpy as np
import matplotlib.patches as patches
import scipy.special as spc

"""
Example of muon's generation with our simulation framework

In the following we are going to evaluate the ionization produced in a given 
medium (Argon, in this case) by muons of a given Energy. This will result in a 
number of clusters and electrons/cluster.
"""
plt.close('all')
# =============================================================================
# We take the Energy from a Poisson distribution centered in ~5 Gev. We can be 
# more precisse in this. (Reference in the main class code)
# =============================================================================

#Inputs -----> from file (FIXME)
n_tracks = 5
E = np.random.poisson(lam = 4000)   #E = 4000 [MeV] (CR)
muons = []                          #list where muon tracks will be stored
dimensions = [15,10,10]
bins = 100
sigma_diff = 0.018
sigma_PSF = 0

#Muons generation
muon = muon_generator(energy = E, geometry = dimensions) #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons) #generate muon's obj stored in muons list

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


"""

muon1, muon2, muon3, muon4, muon5 = muons


plt.figure()

plt.plot(muon1.electron_positions[:,0], muon1.electron_positions[:,1])
plt.plot(muon1.electron_positions[:,0], muon1.electron_positions_diff[:,1], 'o')

plt.plot(muon2.electron_positions[:,0], muon2.electron_positions[:,1])
plt.plot(muon2.electron_positions[:,0], muon2.electron_positions_diff[:,1], 'o')

plt.plot(muon3.electron_positions[:,0], muon3.electron_positions[:,1])
plt.plot(muon3.electron_positions[:,0], muon3.electron_positions_diff[:,1], 'o')

plt.plot(muon4.electron_positions[:,0], muon4.electron_positions[:,1])
plt.plot(muon4.electron_positions[:,0], muon4.electron_positions_diff[:,1], 'o')

plt.plot(muon5.electron_positions[:,0], muon5.electron_positions[:,1])
plt.plot(muon5.electron_positions[:,0], muon5.electron_positions_diff[:,1], 'o')
"""
