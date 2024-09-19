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

gain = 500

magnification = 19.57
diff = 0.20
ath_angle = 130 * deg_to_rad
phi_angle = 47 * deg_to_rad
bins = [int(2688 / rebin),int(2196 / rebin)]
pixel_size = (4.5e-4 * rebin * magnification)                                  #real size * rebinning * magnification 

x_range = [-(pixel_size * bins[0])/2, (pixel_size * bins[0])/2]
y_range = [-(pixel_size * bins[1])/2, (pixel_size * bins[1])/2]

opt_params = LoadOpticalParams(rebin)
#GENERATION ===================================================================

source=Source(radius = 0.35, red_fact = red, rate = 1750)
track_list=[]
"""
exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg", ath_in = ath_angle, phi_in = phi_angle)

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

IMAGENtiff = noise.add_noise(exposition_time, image2d.Hist2D).T
#PLOTS ========================================================================
plt.figure()
plt.title('Optical Gain = {}'.format(gain))
img = plt.imshow(IMAGENtiff, cmap = 'gray', vmax = np.max(IMAGENtiff) * 1.1)
plt.colorbar(img)
"""
#==============================================================================
#==============================================================================
#==============================================================================
# GAIN 1500

gain = 750
track_list=[]

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg", ath_in = ath_angle, phi_in = phi_angle)

diff_handler=Diffusion_handler(sigma_diff = diff, sigma_PSF = 0)
diff_handler.diffuse(track_list)

for i in track_list: i.fill()

noise = Noise(410)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins": bins, "range":[x_range, y_range]}, 
                   QE = float(opt_params['qeff']), GE = float(opt_params['geomeff']), Tp = float(opt_params['T']), 
                   gain = gain)
image2d.track_plot()
image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x(variable_cut = 'x', min_var = -0.648, max_var = 0.378)

IMAGENtiff0 = noise.add_noise(exposition_time, image2d.Hist2D)

plt.figure()
plt.title('Optical Gain = {}'.format(gain))
img = plt.imshow(IMAGENtiff0, cmap = 'RdBu', vmax = np.max(IMAGENtiff0) * 1.1)
plt.colorbar(img)



#==============================================================================
#==============================================================================
#==============================================================================
#GAIN 3000

gain = 3000
track_list=[]

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg", ath_in = ath_angle, phi_in = phi_angle)

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

IMAGENtiff1 = noise.add_noise(exposition_time, image2d.Hist2D)

plt.figure()
plt.title('Optical Gain = {}'.format(gain))
img = plt.imshow(IMAGENtiff1, cmap = 'RdBu', vmax = np.max(IMAGENtiff1) * 1.1)
plt.colorbar(img)

#==============================================================================
#==============================================================================
#==============================================================================
#GAIN 5000

gain = 5000
track_list=[]

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg", ath_in = ath_angle, phi_in = phi_angle)

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

IMAGENtiff2 = noise.add_noise(exposition_time, image2d.Hist2D)

plt.figure()
plt.title('Optical Gain = {}'.format(gain))
img = plt.imshow(IMAGENtiff2, cmap = 'RdBu', vmax = np.max(IMAGENtiff2) * 1.1)
plt.colorbar(img)

#==============================================================================
#==============================================================================
#==============================================================================
#GAIN 10000
"""
gain = 10000
track_list=[]

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg", ath_in = ath_angle, phi_in = phi_angle)

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

plt.figure()
plt.title('Optical Gain = {}'.format(gain))
img = plt.imshow(IMAGENtiff, cmap = 'gray', vmax = np.max(IMAGENtiff) * 1.1)
plt.colorbar(img)
"""

#==============================================================================
#==============================================================================
#==============================================================================
# LOAD DATA FOR COMPARISON

img0 = Image.open('../../Gain anlysis/Data/20240726/c-5254_r-2698_ac-2410_aa-300_fcground_fa1190/12x12/2ms/ss_single_42.tiff')

imgArr = np.array(img0)
#Cutting the blob in the center 

for i in range(79,85):
    for j in range(114,118):
        imgArr[i][j] = np.random.poisson(150) #100

img0 = Image.fromarray(imgArr)

plt.figure()
plt.title('Data from 07-24')
cbar = plt.imshow(img0, cmap = 'RdBu', vmax = np.max(img0) * 1.1)
plt.colorbar(cbar)

#==============================================================================

img1 = Image.open('../tiffs/test/5 ms even higher gain/12x12 pre exposure/ss_single_5.tiff')

imgArr = np.array(img1)

for i in range(95,103):
    for j in range(104,110):
        imgArr[i][j] = 400

img1 = Image.fromarray(imgArr)

plt.figure()
plt.title('Data from 07-23')
cbar = plt.imshow(img1, cmap = 'RdBu', vmax = np.max(img1) * 1.1)
plt.colorbar(cbar)


#==============================================================================
#==============================================================================
#==============================================================================

img = Image.open('../../Gain anlysis/Data/20240726/hvoff/12x12/2ms/ss_single_5.tiff')

plt.figure()
plt.title('Data for bkg')
cbar = plt.imshow(img, cmap = 'RdBu', vmax = np.max(img) * 1.1)
plt.colorbar(cbar)





