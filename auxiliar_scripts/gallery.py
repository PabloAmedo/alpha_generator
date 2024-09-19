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

print('Running...')
#print(os.getcwd())     # --> auxiliar_scritps

path = '../../Gain anlysis/Data/20240719/fr4_acryliclowestpossible/c-3064_r-1928_ac-1800_aa-300_fcground_fa1160/2ms/'
extension = '.tiff'
cmap = 'RdBu'

tiff_files = [f for f in os.listdir(path) if f.endswith('.tiff')]    

for file in tiff_files:
    img = Image.open(path + file)
    
    cmap_change(img, cmap = cmap, title = file)
    #plt.savefig(path + 'cmaps/' + cmap + '/' + file)