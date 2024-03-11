# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 09:49:57 2024

@author: diego
"""

from Alpha_track_simulator import *
import numpy as np
import matplotlib.pyplot as plt
from Clusters_calculation import*

plt.close('all')

#INPUTS
gas1 = 'Argon'
n_tracks = 1000
Emuon = 5000                                                                   #(MeV) (CR)
dimension = [1.4,1.4,1.4]                                                      #Dimension of the chamber
muon_mass = 105.65837                                                          #(MeV/c2)
rho_Ar = 1.784e-3                                                              #(g/cm3)
lines = True                                                                   #True => straight tracks ; False => random tracks

#LOAD DATA (Allison & Cobb paper and electron distribution per cluster)
#Load Allison & Cobb
ACdata = np.loadtxt('TestBoxes.csv', delimiter=';', usecols=1)
bins = len(ACdata)

###############################################################################
###############################################################################

#FUNCTIONS DEFINITION (This could be imported from another .py with aux functions, at least FWHM)
def Resolution_xP(x, n = 1 , P = 1):
    return 96 * n**(-0.46) * (x*P)**(-0.32)

def E_over_I(beta, xP, nu=18, I=188):
    return (6.83*nu*xP)/(I*beta**2)

def FWHM(hval, hbins):
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

###############################################################################
#ARGON=========================================================================

tr_len = []                                                                    #Track lengths storage
most_prob = [] 
E_loss_cm = [] 
e_track_acum = []
clusters = []
e_acum = []
#MUON GENERATION
muons = []
muon = muon_generator(energy = Emuon, geometry = dimension)                    #First generate the muon generator object
muon.produce_muon(n = n_tracks, store = muons, gas = gas1, line = lines)       #generate muon's obj stored in muons list 

#CALCULATIONS
for muon in muons:
    #muon.fill()
    clusters.append(len(muon.n_e_cl))
    tr_len.append(muon.track_length)
    e_track_acum.append(muon.n_electrons)
    E_loss_cm.append((muon.n_electrons * 26.4e-3) / muon.track_length )        #Energy loss / cm (default ~1.5 cm)
    for i in range(len(muon.n_e_cl)):
        e_acum.append(muon.n_e_cl[i])
    

hval, hbins = np.histogram(E_loss_cm, bins = bins)                             #Histograming without showing it

#Calculating FWHM
fwhm = FWHM(hval, hbins)
tr_len = np.mean(tr_len)
total_e = np.mean(e_track_acum)
n_clusters = np.mean(clusters)

#Computing experimental data
E_loss_data, N_cl_data = clusters_cm(Emuon = Emuon)
N_cl_data = 35
# SOME NUMERICAL VALUES =======================================================
print('\n')
print('='*24,'INPUTS','='*24)
print('• Energy:\t\t {}\t\t(MeV)'.format(Emuon))
print('• Momentum:\t\t {:.2f}\t(MeV/c)'.format(np.sqrt(Emuon**2 - muon_mass**2) ))
print('• N tracks:\t\t {}'.format(n_tracks))
print('• Track length:\t {:.2f}\t\t(cm)'.format(muons[0].track_length))
print('• Rho Ar:\t\t {}\t(g/cm3)'.format( rho_Ar))
print('\n')

print('='*15,'RESULTS FROM SIMULATION','='*15)
print('• FWHM:\t\t\t {:.2f}\t\t(keV)'.format(fwhm))
print('• Peak:\t\t\t {:.2f}\t\t(keV)'.format(max(set(E_loss_cm), key=E_loss_cm.count)))
print('• <e-/cm>:\t\t {:.2f}'.format(total_e / tr_len))
print('• <e-/cl>:\t\t {:.2f}'.format(total_e / n_clusters))
print('• <cl/cm>:\t\t {:.2f}'.format(n_clusters / tr_len))
print('• Max e-:\t\t {:.2f}'.format( max(e_acum)))
print('• <dE/dx> sim:\t {:.3f}\t(keV·cm2/g)'.format(np.mean(E_loss_cm) / rho_Ar))
print('\n')

print('='*18,'VALUES FROM DATA','='*18)
print(' '*15,'Data for 3GeV/c pions (-)',' '*15)
print('• FWHM:\t\t\t {:.2f}\t\t(keV)'.format(2.1))                             #From Harris 1972
print('• Peak:\t\t\t {:.2f}\t\t(keV)\t !!!'.format(2.2))                       #From Harris 1972
print('• <e-/cm>:\t\t {:.2f}'.format(3.3 * N_cl_data))
print('• <e-/cl>:\t\t {:.2f}'.format(3.3))                                     #Value from calculations using: e_in_clusters_TEST.py
print('• <cl/cm>:\t\t {:.2f}'.format(N_cl_data))
print('• Max e-:\t\t {:.2f}'.format( 104.0))
print('• <dE/dx> NIST:\t {}\t\t(MeV·cm2/g)'.format(E_loss_data))

###############################################################################
###############################################################################