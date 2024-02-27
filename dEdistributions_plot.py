# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 09:49:57 2024

@author: diego
"""

from Alpha_track_simulator import *
import numpy as np
import matplotlib.pyplot as plt
from Ionization_resolution import FWHM

"""




===============================================================================


THIS IS NOT WORKING PROPERLY !!!!!!!!!!!


===============================================================================




"""

plt.close('all')
# =============================================================================
# We take the Energy from a Poisson distribution centered in ~5 Gev. We can be 
# more precisse in this. (Reference in the main class code)
# =============================================================================

#Inputs -----> from file (FIXME)
n_tracks = 1000
Emuon = 4000   #E = 4000 [MeV] (CR)
muons = []                                                                     #list where muon tracks will be stored
dimensions = ([1.0,1.0,1.0], [1.5,1.5,1.5])#, [2.5,2.5,2.5], [5,5,5], [10,10,10])#, [500,500,500])
### WARNING!! --> Idk why but the previous line ^^^ only works for 3+ dimensions (???)

muon_mass = 105.66 #MeV
samples = 500
###############################################################################

#DATA LOAD -- FROM A&C PAPER
ACdata = np.loadtxt('TestBoxes.csv', delimiter=';', usecols=1)
xaxis_data = np.linspace(0, 10, len(ACdata))
bins = len(ACdata)
"""
We don't really care about the diffusion and noise, so we get rid out of them. 
We are just interested in the total number of e- generated and the energy loss.
"""
E_peak = []
widths = []
dE_distribution = []

#Now we want to separate the traze in a given number of samples and calculate 
#the amount of energy in each one.
for dimension in dimensions:
    #MUON GENERATION
    muons=[]
    muon = muon_generator(energy = Emuon, geometry = dimension)                #First generate the muon object
    muon.produce_muon(n = n_tracks, store = muons, line=True)                             #generate muon's obj stored in muons list
    
    for muon in muons:
        dE_distribution = []
        tr_len = muon.track_length
        n_e_cl = np.array(muon.n_e_cl)
        n_clusters = len(muon.n_e_cl)
        n_electrons = muon.n_electrons
        sep =  45  
        suma = 0                                                            #this is causing problems in < 1 cm tracks
        for i in range(0, n_clusters, sep):  
            for j in range(i, i + sep):
                if j < n_clusters:
                    suma = suma + muon.n_e_cl[j]
                else:
                    break
                #except IndexError as err:
                 #   suma = np.mean(dE_distribution) / 26.4e-3
                    #break
            dE_distribution.append(suma * 26.4e-3) #keV
    
    hval, hbins, _ = plt.hist(dE_distribution, bins=bins)
    plt.figure()
    hval, hbins, _ = plt.hist(dE_distribution, bins=bins)
    peak = np.argmax(hval)
    #fwhm = FWHM(hval, hbins)
    plt.title('Track legth: {:.2f} (cm) // ({}) tracks'.format(muon.track_length, n_tracks))
    plt.plot(xaxis_data, ACdata/np.sum(ACdata), 'o', color = 'r', label = 'exp data [A&C]')
    #plt.plot([],[], ' ',  label = 'FWHM: {}'.format(fwhm))
    plt.plot([],[], ' ',  label = 'Peak: {} (keV)'.format(hbins[peak]))
    plt.xlabel('Energy loss (keV)')
    plt.ylabel('Counts')
    plt.legend()





