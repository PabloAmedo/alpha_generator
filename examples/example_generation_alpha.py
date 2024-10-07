# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 11:16:45 2023

@author: Pablo
"""
import os
os.chdir('../')

from Alpha_track_simulator import*
import sys
import time

sys.path.append('/..')

st = time.time()
# This is an example of how to generate alpha tracks using this generator
print('Running...')
# Create an argon object from the gas class
argon = Gas(gas = 'ArCF4')

# Create a source and a list to store the tracks

source = Source(radius = 0.35, red_fact = 100) #Radius in cm

track_list = []

# Set some parameters for the alpha track production
n_tracks = 1 ; ath_angle = 0

exposition_time = n_tracks / source.rate #In seconds

# Generate the alpha tracks using the Bragg profile
source.produce_alpha(n = n_tracks, store = track_list, ath_in = None, ionization_profile = "Bragg")


#Create a difussion handler and change the diffusion of the tracks (sigma in cm)
diff_handler = Diffusion_handler(sigma_PSF = 0)
# Update the diffusion for each track
diff_handler.diffuse(track_list, 180)
"""
# Generate the electrons alongside the tracks. This will generate the original positions of
# the electrons and the diffused positions as well
for i in track_list: i.fill(diff = True)


#Create a noise object with a given dark noise
noise = Noise(50)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins":100})
# #Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x()
"""
end = time.time()

print('Elapsed time:\t', end - st)









