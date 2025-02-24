# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:34:44 2024

@author: diego
"""

import numpy as np
import os
import time
from PIL import Image
from PIL.ExifTags import TAGS
from functions import *
import matplotlib.pyplot as plt

os.chdir('../')
from Alpha_track_simulator import *

st = time.time()

#INPUTS =======================================================================
n_tracks    =   1750//5
red         =   1000
native_bin  =   4
P           =   2.5

diff        =   0.22 

bins        =   [int(4096 / native_bin),int(2304 / native_bin)]                #Change to binning --> this is the size for a 12x12 reb
pixel_size  =   (4.6e-4 * native_bin * 19.48 )                           #real size * rebinning * magnification 
QE          =   0.7
GE          =   1.6511e-4
T           =   1

path        =   '../../NOV24 - EXPERIMENTAL CAMPAIGN/Camera/20122024/Both/200ms/S1/'

x_range     =   [-(pixel_size * bins[0])/2, (pixel_size * bins[0])/2]
y_range     =   [-(pixel_size * bins[1])/2, (pixel_size * bins[1])/2]

save        =   0
#Cuts in x axis (pxs)
cut1 = 350
cut2 = 650 

print(os.getcwd())

#SIMULATE EXPERIMENTAL CONDITIONS  ============================================
source = Source(radius = 0.35, red_fact = red, rate = 1750, P = P)
track_list = []

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg")

diff_handler=Diffusion_handler(sigma_diff = diff, sigma_PSF = 0)
diff_handler.diffuse(track_list)

for i in track_list: i.fill(diff = True)

noise = Noise(0.06)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins": bins, "range":[x_range, y_range]}, 
                   QE = QE, GE = GE, Tp = T, 
                   gain = 1)

hist, ecut = image2d.plot_x(variable_cut = 'x')#, min_var = -0.648, max_var = 0.378)
primaries = image2d.Hist2D_e
primaries_profile = np.sum(primaries, axis = 0) * red

recorded = image2d.Hist2D * red

#LOAD DATA + GAIN =============================================================
data_prof, gains, gains2, profs = GetGain(path, 
                           primaries_profile = primaries_profile,
                           cuts = ([138, 458], [cut1, cut2]),
                           QE = QE, GE = GE, T = T,
                           show = True)

#PLOTS AND RESULTS ============================================================
plt.figure()
plt.hist(gains, bins = 50, histtype = 'step', label = 'Peaks')
plt.hist(gains2, bins = 25, histtype = 'step', label = 'Integral')
plt.xlabel('Optical Gain (phe/e)', fontsize = 15)
plt.legend()

x = np.linspace(0, len(data_prof), len(data_prof))

sec = image2d.Hist2D * red / QE / GE / T 
sec_profile = np.sum(sec, axis = 0)
meanprofs = np.mean(profs, axis = 0)

# Plot profile comp
"""
plt.figure()
plt.plot(x, profs / max(profs),  label = 'Data')
plt.plot(x - 7, primaries_profile / max(primaries_profile), label = 'Sim')
plt.xlabel('px', fontsize = 15)
plt.ylabel('phe', fontsize = 15)
plt.xlim([cut1, cut2])
plt.legend()
"""

plt.figure()
plt.plot(meanprofs/max(meanprofs[400:550])) #200:800
plt.plot(x-20 + 12, primaries_profile / max(primaries_profile))
plt.plot(x[340:350] ,meanprofs[355:365]/max(meanprofs[400:550])) #200:800


print('='*50)
print('Averaged over:\t\t\t{:.3f} files'.format(len(gains)))
print('Aerage peaks:\t\t\t{:.3f} (phe/e)'.format(np.mean(gains)))
print('Aerage integral:\t\t{:.3f} (phe/e)'.format(np.mean(gains2)))

print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(max(gains)))
print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(min(gains)))

end = time.time()

print('\nElapsed time:\t{} s'.format(end - st))


