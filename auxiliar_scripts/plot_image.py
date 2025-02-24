# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:41:28 2024

@author: usuario
"""

import numpy as np
import matplotlib.pyplot as plt
from optical_gain_tools import *
from PIL import Image

folder2 = '../tiffs/July24/highgain_compare/10ms/'
path2 = folder2 + 'ss_single_1'

folder = '../tiffs/test/bin 12x12 5ms/'
path = folder + 'ss_single_13'

ext = '.tiff'

base = 108
bkg1 = 0.35
bkg2 = 0.35

data1__ = Image.open(path + ext)
data2__ = Image.open(path2 + ext)

data1 = BaselineSubstraction(data1__, baseline = 412)
data2 = BaselineSubstraction(data2__, baseline = 108)

fig, ax = plt.subplots(1, 2, figsize = (18,18))

fig.suptitle('Acrylic only', fontsize = 25)

ax[0].imshow(data1, cmap='jet')
ax[0].set_title('July 2023 (5 ms)', fontsize = 20)
ax[0].set_xlabel('X (px)', fontsize = 18)
ax[0].tick_params(axis='both', which='major', labelsize = 16)

ax[1].imshow(data2, cmap='jet')
ax[1].set_title('July 2024 (10 ms)', fontsize = 20)
ax[1].set_xlabel('X (px)', fontsize = 18)
ax[1].tick_params(axis='both', which='major', labelsize = 16)

fig.tight_layout()

plt.show()





