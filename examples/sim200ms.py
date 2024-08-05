# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:12:12 2024

@author: diego
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:09:28 2024

@author: diego
"""

import os
os.chdir('../')

from PIL import Image
import tifffile
from Alpha_track_simulator import*
from optical_gain import*
import sys
import time

print('Running...')
sys.path.append('/..')

start = time.time()
plt.close('all')

#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)                                                           #THIS IS NOT BEING USED AT ALL
#Create a source and a list to store

source=Source(radius = 0.35, red_fact = 1)
track_list=[]

#Creat some alpha tracks and set some parameters
n_tracks = 100 ; ath_angle = 0

exposition_time = n_tracks/source.rate #In seconds
source.produce_alpha(n = n_tracks, store = track_list, ath_in = None, ionization_profile = "Bragg")

#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler(sigma_diff=0.025,sigma_PSF=0)
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()

#Create a noise object with a given dark noise
noise = Noise(50)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins":100})          #Calcualted to match 12x12 rebinning px size (4.5*12 (um))
# #Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x(variable_cut = 'x', min_var = -0.648, max_var = 0.378)
os.chdir('alpha_generator/')
plt.figure()
plt.imshow(image2d.Hist2D)

np.savetxt('../sim200ms.txt', image2d.Hist2D)