# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 14:06:04 2024

@author: diego
"""
print('Running...')

import numpy as np

def GeometricalEfficiency(N, f, distance):
    
    magnification = f / (f - distance)
    diameter = f / N
    e_geo = np.pi * (diameter / 2)**2 / (4 * np.pi * (f * (1 - 1/magnification))**2)
    return magnification, e_geo

m, eg = GeometricalEfficiency(0.95, 2.5, 51.9)
print(eg)