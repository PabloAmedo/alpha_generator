# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:54:23 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from optical_gain_tools import *
import os
from PIL import Image

print('Running...')

path = '../../Gain anlysis/data_aux/c-5254_r-2698_ac-2410_aa-300_fcground_fa1190/12x12/100ms/ss_single_1.tiff'
img = Image.open(path)

image = BaselineSubstraction(img, baseline=110)

xproj = np.sum(image, axis = 0)

