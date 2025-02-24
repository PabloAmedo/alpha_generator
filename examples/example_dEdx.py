# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 16:21:08 2024

@author: diego
"""

import os
os.chdir('../')
from Alpha_track_simulator import *
import time
st = time.time()

# INPUTS ======================================================================
n_tracks = 10000
dimensions = [1.5, 1.5, 1.5] #cm
gas = 'ArCH4_93-7'
P = 10 #bar
mass = 139.57 #MeV/c2 --> pi-
energy = 3500
#e_cut = 100000
# =============================================================================
print(os.getcwd())
os.chdir('alpha_generator/')
beam =  muon_generator(energy = energy, geometry = dimensions, gas= gas, mass = mass, pressure = P)
ne = []
for i in range(n_tracks):
    track_list = []
    beam.produce_muon(n = 1, store = track_list, line=True)
    ne.append(track_list[0].n_electrons)

ne = np.array(ne)
dE = ne * 26.7e-3 #keV

data_Harris = np.loadtxt('data/Harris72.csv', delimiter = ';')

hval, hbin = np.histogram(dE, bins = 75)

p80 = np.percentile(dE, 87.5)
hval, hbin, patches = plt.hist(dE, bins=50, alpha=0.7)

# Resaltar los valores mayores al percentil 80
for i in range(len(hbin) - 1):
    if hbin[i] >= p80:
        patches[i].set_facecolor('gray')  # Cambia al color deseado

# Mostrar la gráfica
plt.axvline(hbin[hbin>p80][0], color='red', linestyle='solid', label=f'80th percentile ({p80:.2f})')
plt.show()
"""
plt.figure()
#plt.title('Ecut = {:.2f} (MeV)'.format(e_cut * 26.7e-6))
#plt.plot(data_Harris[:,0], (data_Harris[:,1]/max(data_Harris[:,1])) * max(hval), marker = 'x', linestyle = '-.', label = 'Data', alpha=0.7, 
#         color='blue', linewidth=1.5)
plt.hist(dE*1.07, bins = 50, histtype = 'stepfilled', label = 'Simulated', alpha=0.5, color='orange', edgecolor='black')

#plt.legend()
"""
end = time.time()
print('Elapsed time:\t{:.2f}'.format(end - st))
    

