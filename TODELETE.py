# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:29 2024

@author: diego
"""

from Alpha_track_simulator import *
import numpy as np
import matplotlib.pyplot as plt


#Avoid warning in console
import warnings
from scipy.optimize import OptimizeWarning
warnings.filterwarnings("ignore", message="Covariance of the parameters could not be estimated", category=OptimizeWarning)


#INPUTS
gas = 'Argon'
n_tracks = 1000
Emuon = 3000                                    #(MeV) (CR)
mass_muon = 139.57                              #(MeV/c2)  PI
p = np.sqrt((Emuon + mass_muon)**2 - mass_muon**2)
dimensions = ([1.45,1.45,1.45])
muon_mass = 139.57                              #(MeV) for beta calculation
rho_Ar = 1.784                                  #rho Ar/CH4 (93/7)
W_Ar = 26.4e-3                                  #keV


#Load Allison & Cobb
ACdata = np.loadtxt('TestBoxes.csv', delimiter=';', usecols=1)
bins = len(ACdata) * 3                         #I think this is around 64 * 2 bins
x_AC_data = np.linspace(0, 10, len(ACdata))
#Generate lcp
muons = []
muon = muon_generator(energy = Emuon, geometry = dimensions)                   #First generate the muon generator object
muon.produce_muon(n = n_tracks, store = muons, gas = gas, line = True, e_cut=775)

E_loss_total = []
E_loss_cm = []

for track in muons:
    E_loss = track.n_electrons * W_Ar
    length = track.track_length
    
    E_loss_total.append(E_loss)
    E_loss_cm.append(E_loss / length)

E_loss_total = np.array(E_loss_total)
E_loss_cm = np.array(E_loss_cm)


"""
plt.figure()
plt.title('Total E loss')
hval, _,_ =plt.hist(E_loss_total / ( rho_Ar), bins = bins)
plt.plot(x_AC_data, ACdata/(max(ACdata)/max(hval)), 'kx')
plt.plot(3.95, max(hval), 'ro')
"""


plt.figure()
plt.title('E loss per unit length')
hval, hbins,_ =plt.hist(E_loss_cm   , bins = bins, label = 'simulated')
plt.plot(x_AC_data, ACdata/(max(ACdata)/max(hval)), 'kx', label = 'Harris72 data')

plt.legend()

    