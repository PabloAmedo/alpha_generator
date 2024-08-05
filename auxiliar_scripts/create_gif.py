# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:46:21 2024

@author: diego
"""

import imageio
import re
import numpy as np
import matplotlib.pyplot as plt
import os

print('Running...')

os.chdir('../')

rebin = '12x12_test'                                                             #based on experimental settings
exposition_time = 1 #s                                                       #based on experimental settings
path_optData = 'alpha_generator/setup_params/'                                 #optical parameters stored in a .txt file
#folder = 'tiffs/opt gain/optical_gain_10-07-24/fr4 on 2/1 s/'#+ str(exposition_time) +' s/'          #images location --> all stored in one folder
folder = 'tiffs/July24/new set data/highgain-2structures/12x12/100ms/cmap2/'                           

tiff_files = [f for f in os.listdir(folder) if f.endswith('.tif')]
def extract_number(file_name):
    match = re.search(r'(\d+)', file_name)
    return int(match.group(1)) if match else 0

sorted_tiff_files = sorted(tiff_files, key = extract_number)

images = []
for tiff in sorted_tiff_files:
    images.append(imageio.imread(folder + tiff))
imageio.mimsave('100ms_test_5k.gif', images, duration = 0.25)