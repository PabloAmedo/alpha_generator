# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:32:36 2024

@author: diego
"""

print('Running...')

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import os

print(os.getcwd())
"""
offset = 100
data_img = Image.open('../tiffs/July/highgain_compare/1000ms/ss_single_1.tiff')
#data_img = Image.open('../../../../CCD CAMERA/Noise measurements/noise_08-05-2024/lights off/10ms/12x12.tiff')
data = np.array(data_img)

data[data < offset] = 0
data[data >= offset] -= offset
print(np.sum(data))
"""
A = 500 #Bq
W = 25.7e-6 #MeV
E = 5.5 #MeV

geomeff = 0.0001967737393285824
qeff = 0.7
T = 1

def ExpectedPhotons(time_exp):
    return time_exp * geomeff * qeff * T * A * E / W

times = np.array([ 0.05, 0.100, 0.200, 0.700, 1])
gains = np.array([ , , 79.696, , , ])

expected_photons = ExpectedPhotons(times)

plt.figure()
plt.title('Gain vs exposure time', fontsize = 25)
plt.plot(times, gains, 'k-', lw = 2, label = 'Gain from script')
plt.xlabel('Exposure time (s)', fontsize = 20)
plt.ylabel('Gain (ph/primary e)', fontsize = 20)

