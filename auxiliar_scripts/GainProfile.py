# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:04:16 2024

@author: diego
"""
import os
print('Running...')
os.chdir('../')

from Alpha_track_simulator import*
import tifffile
from optical_gain import*
import pandas as pd
from scipy.optimize import curve_fit
from tools import *
import time

initial = time.time()
os.chdir('alpha_generator/') #!!! THIS SHOULD BE A GLOBAL PATH IF YOU DOWNLOWAD FROM GITHUB!!!
#INPUTS =======================================================================
simulation      = 0         #Select if you want to run a new simulation (set to 1) or use stored data (set to 0)
f               = 500       #Bq
#exposition_time = 0.5      #s
offset          = 100       #Camera ADU's baseline/offset
bins            = 100       
bkg_cut         = 0.05      #percentage of the peak where we cut to consider bkg --> THIS COULD BE OPTIMIZED!!
rebinning       = '12x12'
last_structure  = 'test' #'FR4on'
#other stuff needed
n               = 100 #int(f * exposition_time) #number of tracks generated based on exp time
exposition_time = n / f
red_fact        = 10
alphas          = []                       #list to store the tracks

#LOAD DATA ====================================================================
#Load optical parameters
opt_params      = pd.read_csv('setup_params/datos' + rebinning + '_' + last_structure + '.csv')

#Load camera image
set_path        = 'tiffs/'
run             = 'July24/highgain_compare/'
exp_time        = '200ms/'
file            = 'ss_single_1.tiff'
image_tiff      = Image.open(set_path +  run + exp_time+file)
data_image      = np.array(image_tiff)

#GENERATION / LOAD ============================================================
if simulation == 1: #GENERATION
    #Source definition + alpha production
    source = Source(radius = 0.35, red_fact = red_fact)
    source.produce_alpha(n = n, store = alphas, ionization_profile = "Bragg")
    #Diffusion
    diff_handler = Diffusion_handler(sigma_diff = 0.25, sigma_PSF = 0)
    diff_handler.diffuse(alphas)
    #Fill tracks with electrons
    for i in alphas: i.fill()
    #Noise
    noise = Noise(50) #Check this value
    #Image formation
    image2d = Image_2D(track_list = alphas, hist_args = {"bins":bins})         #This defines the image2d object
    image2d.plot_hist(noise_object = noise, exposition_time = exposition_time) #Plots the simulated image
    image2d_noise = noise.add_noise(exposition_time, image2d.Hist2D)
    image2d.plot_x()
else: #LOAD
    data_sim = np.loadtxt('data/simulated_alphas/Simulated1s')


#data_image = HotPixel12x12(data_image)
data_substracted = BkgSubstraction(data_image, show = True)

x = np.linspace (1, len(data_substracted), len(data_substracted)) - np.argmax(data_substracted)

#SHOW AGREEMENT ===============================================================
plt.figure()
plt.title('Simulated vs Data comparison', fontsize = 25)
if simulation ==1:
    hval, hbins, _ = plt.hist(image2d.electron_cut, bins = bins, label = 'Simulated', histtype = 'step')
else:
    hval, hbins, _ = plt.hist(data_sim, bins = bins, label = 'Simulated', histtype = 'step')
plt.plot(x / float(opt_params['cal_fab']), data_substracted / max(data_substracted) * max(hval), label = 'Data {:.3f} (s)'.format(exposition_time))
plt.xlabel('X (cm)', fontsize = 20)
plt.ylabel('n photons', fontsize = 20)
plt.vlines(-4.594315, 0, max(hval) + max(hval)*0.15, 'k', label = 'Alpha range')
plt.vlines(4.594315, 0, max(hval) + max(hval)*0.15, 'k')
plt.ylim([0, max(hval) + max(hval)*0.15])
plt.legend(fontsize = 18)
plt.grid()

final = time.time()
print('Time elapsed:\t{:.2f} (s)'.format(final-initial))
#GAIN =========================================================================
"""
#gain_img = Image.fromarray(data_image)
#gain_img = np.array(data_image)

gain = gainFunc(data = data_image, 
                qeff = float(opt_params['qeff']), 
                geomeff = float(opt_params['geomeff']), 
                T = 1, 
                exp_time = exposition_time,
                base = offset)
"""
