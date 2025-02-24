# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:58:48 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from general_tools import *
import os
from PIL import Image
os.chdir('../../SingleTrackAnalysis-Liv23/rebin')
from rebinning_tools import *
print(os.getcwd())
os.chdir('../../alpha_generator/auxiliar_scripts')

print('Running...')
print(os.getcwd())
path = '../../../../NOV24 - EXPERIMENTAL CAMPAIGN/Camera/20122024/Both/10ms/'
extension = '.tiff'
cmap = 'RdBu'


tiff_files = [f for f in os.listdir(path) if f.endswith('.tif')]    

for file in tiff_files:
    img_ = Image.open(path + file)
    
    imgar = np.array(img_, dtype = np.uint8)
    imgreb = rebin(imgar[:,:-1], 3)
    
    
    cmap_change(imgreb, cmap = cmap, title = file, pct_max = 1)
    plt.savefig(path + 'cmaps/' + cmap + file)