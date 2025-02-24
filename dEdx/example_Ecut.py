# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:52:50 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
"""
Example to visualize the impact of the Ecut in the dEdx simulations.
All the values have been obtained through previous simulations for the ND-GAr 
settings:
    - 5m
    - 10 bar
    - ArCF4 99/1
    - muons in MIP region (E ~ 250 MeV)
"""

Ecut = np.array([0.267, 2.67, 13.35, 26.7, 128])                #MeV
dEdx_avg = np.array([8.961, 10.669, 11.859, 12.45, 13.178])     #MeV/5m
dEdx_trunc = np.array([8.96, 10.45, 11.03, 11.1, 11.28])        #MeV/5m accepting 80%

plt.figure()
plt.title('E primary = 242 (MeV)')
plt.plot(Ecut, dEdx_avg, 'b', label = 'Average')
plt.plot(Ecut, dEdx_trunc, 'r', label = 'Truncated')
plt.xlabel('E-cut (MeV)', fontsize = 15)
plt.ylabel('<dE> (MeV/5m)', fontsize = 15)
plt.legend()