# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 09:49:57 2024

@author: diego
"""

from Alpha_track_simulator import *
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
# =============================================================================
#FIXME
#check for missnames
# =============================================================================

#INPUTS
gas1 = 'Argon'
gas2 = 'Xenon'
gases = (gas1, gas2)
n_tracks = 50000
Emuon = 2285                        #(MeV) (CR)
dimensions = ([0.3,0.3,0.3], [1.4,1.4,1.4])    #Different geometries for the track lengths (format is weird, I know)
bins = 64                           #bins for the histogram
muon_mass = 105.66                  #(MeV) for beta calculation
rho_Ar = 1.784
norm_index = 1

#LOAD DATA (Allison & Cobb paper and electron distribution per cluster)
#Load Allison & Cobb
ACdata = np.loadtxt('TestBoxes.csv', delimiter=';', usecols=1)
bins = len(ACdata)
###############################################################################

#FUNCTIONS DEFINITION (This could be imported from another .py with aux functions, at least FWHM)
def Resolution_xP(x, n = 1 , P = 1):
    return 96 * n**(-0.46) * (x*P)**(-0.32)

def E_over_I(beta, xP, nu=18, I=188):   #Not sure at all about this
    return (6.83*nu*xP)/(I*beta**2)

def FWHM(histogram, bins):
    # Encuentra la posición del valor máximo del histograma
    max_index = np.argmax(hval)
    max_value = hval[max_index]
    max_position = (hbins[max_index] + hbins[max_index + 1]) / 2.0
    # Encuentra la mitad del valor máximo
    half_max_value = max_value / 2.0
    # Encuentra los índices más cercanos al lado izquierdo y derecho de la mitad del valor máximo
    left_index = np.argmin(np.abs(hval[:max_index] - half_max_value))
    right_index = np.argmin(np.abs(hval[max_index:] - half_max_value)) + max_index
    # Calcula FWHM
    fwhm = hbins[right_index] - hbins[left_index]
    
    return fwhm
    
#ARGON=========================================================================
tr_len_Ar = []                         #Track lengths storage
widths_Ar = []
most_prob = [] 
#IONIZATION RESOLUTION CALCULATIONS
for dimension in dimensions:
    muons = []
    muon = muon_generator(energy = Emuon, geometry = dimension)                #First generate the muon generator object
    muon.produce_muon(n = n_tracks, store = muons, gas = gas1, line = True)                 #generate muon's obj stored in muons list
    beta = np.sqrt(1 - 1 / (Emuon / muon_mass + 1)**2)                         #this is only needed if we want to use the dimensionless detector length
    E_loss_cm = []                                                             #basically de dEdx
    for muon in muons:
        #muon.fill()
        length = ( muon.track_length)
        n_electrons = muon.n_electrons
        E_loss_cm.append((n_electrons * 26.4e-3) / length )   #eV
    
    """
    We can hide the next histograms if we do not want to plot them by calling 
    np.histogram instead plt.hist (it does exactly the same)
    """
    #E_loss_cm = (np.array(E_loss_cm) / length)
    plt.figure()    #This will show a plot for every dimension we have -------> Must be commented
    hval, hbins,_ = plt.hist(E_loss_cm, bins = bins)
    #x_AC_data = np.linspace(min(hbins), max(hbins), len(ACdata))
    x_AC_data = np.linspace(0, 10, len(ACdata)) + (min(hbins))
    plt.plot(x_AC_data, ACdata/(max(ACdata)/max(hval)), 'o', color = 'r', label = 'exp data [A&C]')
    plt.xlabel('dE / dx (keV / cm)')
    plt.ylabel('Counts')
    plt.title('Track legth: {:.2f} in {}'.format(muon.track_length, gas1))

    hval, hbins = np.histogram(E_loss_cm, bins = bins)                         #Histograming without showing it
    #Calculating FWHM
    fwhm = FWHM(hval, hbins)
    widths_Ar.append(fwhm)
    tr_len_Ar.append(length)
    print('\n')
    print('FWHM:\t', fwhm)
    print('Peak:\t', max(set(E_loss_cm), key=E_loss_cm.count))

tr_len_Ar = np.array(tr_len_Ar)
widths_Ar = np.array(widths_Ar)
widths_norm_Ar = widths_Ar / widths_Ar[norm_index]
E_I_Ar = E_over_I(beta = beta, xP = tr_len_Ar)

#XENON=========================================================================
"""
#This was just used for the verification with Xe too ------ Uncomment if needed

tr_len_Xe = []                         #Track lengths storage
widths_Xe = [] 
#IONIZATION RESOLUTION CALCULATIONS
for dimension in dimensions:
    muons = []
    muon = muon_generator(energy = Emuon, geometry = dimension)                #First generate the muon generator object
    muon.produce_muon(n = n_tracks, store = muons, gas = gas2)                             #generate muon's obj stored in muons list
    beta = np.sqrt(1 - 1 / (Emuon / muon_mass + 1)**2)                         #this is only needed if we want to use the dimensionless detector length
    E_loss_cm = []                                                             #basically de dEdx
    for muon in muons:
        #muon.fill()
        length = ( muon.track_length)
        n_electrons = muon.n_electrons
        E_loss_cm.append((n_electrons * 26.4) / length)   #eV
        
# =============================================================================
#     plt.figure()    #This will show a plot for every dimension we have
#     hval, hbins,_ = plt.hist(E_loss_cm, bins = bins)
#     plt.xlabel('Energy loss per cm (eV/cm)')
#     plt.ylabel('Counts')
#     plt.title('Track legth: {:.2f} in {}'.format(muon.track_length, gas2))
# =============================================================================
    hval, hbins = np.histogram(E_loss_cm, bins = bins)
    #Calculating FWHM
    fwhm = FWHM(hval, hbins)
    widths_Xe.append(fwhm)
    tr_len_Xe.append(length)
    

tr_len_Xe = np.array(tr_len_Xe)
widths_Xe = np.array(widths_Xe)
widths_norm_Xe = widths_Xe / widths_Xe[norm_index]  
E_I_Xe = E_over_I(beta = beta, xP = tr_len_Xe)
"""
#Show the Figure (ENERGY DISTRIBUTION RESOLUTION (?) )=========================
plt.figure()
plt.plot(E_I_Ar, widths_norm_Ar *100 , 'o', color = 'g', label = 'Ar')
#plt.plot(E_I_Xe, widths_norm_Xe *100 , 'x', color = 'r', label = 'Xe')        #Xenon related
plt.plot(E_I_Ar, Resolution_xP(x = E_I_Ar) , '-', color = 'b')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Detector length (E/I adim)')
plt.ylabel('FWHM (%)')
plt.title('Ionization resolution. (tracks={})'.format(n_tracks))
plt.ylim([1,250])
plt.grid()
plt.legend(loc = 'lower left')
