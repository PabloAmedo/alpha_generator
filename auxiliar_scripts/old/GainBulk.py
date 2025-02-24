# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 13:46:30 2023

@author: diego
"""

import os
print('Running...')
os.chdir('../')

from Alpha_track_simulator import*
import tifffile
from optical_gain import*
from PIL import Image
#Inputs

rebin = '12x12_test'                                                             #based on experimental settings
exposition_time = 100 #s                                                       #based on experimental settings
path_optData = 'alpha_generator/setup_params/'                                 #optical parameters stored in a .txt file
#folder = 'tiffs/opt gain/optical_gain_10-07-24/fr4 on 2/1 s/'#+ str(exposition_time) +' s/'          #images location --> all stored in one folder
folder = '../Gain analysis/Data/July11th/driftscan/c-1800_r-900_acground_aa2000/100ms/'                           

datos = pd.read_csv( path_optData + 'datos' + rebin + '.csv' )
os.chdir('alpha_generator/')

tiff_files = [f for f in os.listdir(folder) if f.endswith('.tiff')]            #select speciffic format for images (.tiff)
print('Exp time:\t{:.2f} (s)'.format(exposition_time))
print('Rebinning:\t{}'.format(rebin))
print('-'*50)
gain_list = []
for file in tiff_files:                                                        #iteration over each file
    
    cameraImage = data_image(folder + file)
    print('File: ', file)
    
    gain_list.append(cameraImage.gain(qeff = float(datos['qeff']),
                     geomeff = float(datos['geomeff']),
                     T = float(datos['T']),  
                     exp_time = exposition_time))
    print('-'*50)

meanagain = np.mean(gain_list)
print('Avg Gain:\t', meanagain)