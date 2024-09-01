# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 10:43:20 2024

@author: usuario
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from optical_gain_tools import *
from general_tools import *

plt.close('all')

path = '../../Gain analysis/Data/selected tracks/two tracks/ss_single_4.tiff'
pct = 1.1
img_ = Image.open(path)

img = BaselineSubstraction(img_, baseline = 111)
img_array = np.array(img)
row, col = np.unravel_index(np.argmax(img_array), img_array.shape)

plt.figure()
plt.imshow(img, cmap = 'seismic', vmax = np.max(img) * 1.1)
plt.plot(col, row, 'kx', lw  = 3)


cmap_change(img, 'seismic', pct_max = pct)
cmap_change(img, 'bwr', pct_max = pct)
cmap_change(img, 'jet', pct_max = pct)
cmap_change(img, 'RdBu', pct_max = pct)
cmap_change(img, 'Spectral', pct_max = pct)
cmap_change(img, 'brg', pct_max = pct)
cmap_change(img, 'Blues_r', pct_max = pct)


