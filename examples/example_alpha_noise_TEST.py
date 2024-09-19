# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 13:44:51 2024

@author: diego
"""

import os
os.chdir('../')

from PIL import Image
from Alpha_track_simulator import*
from optical_gain import*

os.chdir('alpha_generator/auxiliar_scripts/')

from optical_gain_tools import *
from general_tools import * 

print('Running...')

#INPUTS =======================================================================

n_tracks = 1
red = 10
rebin = 12
gain = 3500

diff = 0.20 / 3
ath_angle = 220 * deg_to_rad
phi_angle = -45 * deg_to_rad
bins = [int(2688 / rebin),int(2196 / rebin)]                                                                #Change to binning --> this is the size for a 12x12 reb
pixel_size = (4.5e-4 * rebin * 19.57)                                              #real size * rebinning * magnification 

x_range = [-(pixel_size * bins[0])/2, (pixel_size * bins[0])/2]
y_range = [-(pixel_size * bins[1])/2, (pixel_size * bins[1])/2]

opt_params = LoadOpticalParams(rebin)
#GENERATION ===================================================================

source=Source(radius = 0.35, red_fact = red, rate = 1750)
track_list=[]

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg")

diff_handler=Diffusion_handler(sigma_diff = diff, sigma_PSF = 0)
diff_handler.diffuse(track_list)

for i in track_list: i.fill()

noise = Noise(110)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins": bins, "range":[x_range, y_range]}, 
                   QE = float(opt_params['qeff']), GE = float(opt_params['geomeff']), Tp = float(opt_params['T']), 
                   gain = gain)
image2d.track_plot()
image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x(variable_cut = 'x', min_var = -0.648, max_var = 0.378)

IMAGENtiff = noise.add_noise(exposition_time, image2d.Hist2D)
#PLOTS ========================================================================
plt.figure()
plt.title('Optical Gain = {}'.format(gain))
img = plt.imshow(IMAGENtiff, cmap = 'gray', vmax = np.max(IMAGENtiff) * 1.1)
plt.colorbar(img)



