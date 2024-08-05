# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:33:28 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt

time12acr  = np.array([0.1, 0.5, 1, 5, 10])
time8acr  = np.array([0.1, 0.5, 1, 5])
time8fr  = np.array([0.1, 0.3, 0.5, 1])
time12acr_09 = np.array([5, 10, 30])

gain12acr = np.array([1327.1, 662.7, 586.98, 526.8, 512.8]) 
gain8acr = np.array([2063, 1314.74, 665.13, 543.61])
gain8fr = np.array([6401.31, 5867.6, 5439.75, 5425.47])
gain12acr_09 = np.array([273.36, 265.22, 250.18])

plt.figure()
plt.title('Gain vs exposure time and binning', fontsize = 25)
plt.plot(time8acr, gain8acr, 'o-', label = '8x8 - acrylic only (Wed)')
plt.plot(time12acr, gain12acr, 'o-', label = '12x12 - acrylic only (Wed)')
plt.plot(time8fr, gain8fr, 'o-', label = '8x8 - FR4 (Wed)')
plt.plot(time12acr_09, gain12acr_09, 'o-', label = '12x12 - acrylic only (Tue)')
plt.plot(0.1, 1434, 'rv', label = 'Calculation by Jacobo')
plt.xlabel('exposure time (s)', fontsize = 20)
plt.ylabel('gain (ph/e)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize = 16)
plt.legend(fontsize = 16)