# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:10:49 2024

@author: diego
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

print('Running...')

os.chdir('../')
def ClusterParametrizationAr(n):
    return 0.216/n**2 

def ClusterParametrizationCH4(n):
    return 0.119/n**2


e_cut = 10000

#Load Ar/CH4(90/10) data (from HEED presentation)
data_ArCF4_9901_le = np.loadtxt('data/clusters_distributions/ArCF4-99-1-10bar_0.1x0.1x0.1mmCell_0.2GeVmuons_normalized.txt' , skiprows = 1)
n_el_ArCF4_991_le, P_el_ArCF4_991_le, _ = np.split(data_ArCF4_9901_le, 3, axis = 1)

#Load Ar/CF4(99/1) data (from HEED - Carlos)
data_ArCF4_9901_he = np.loadtxt('data/clusters_distributions/ArCF4-99-1-10bar_0.1x0.1x0.1mmCell_8GeVmuons_normalized.txt' , skiprows = 1)
n_el_ArCF4_991_he, P_el_ArCF4_991_he, _ = np.split(data_ArCF4_9901_he, 3, axis = 1)
"""
#Load Pure Ar
data_Ar = np.loadtxt('data/clusters_distributions/PureArgon-10bar_0.1x0.1x0.1mmCell_pi2.5GeV_normalized.txt' , skiprows = 1)
n_el_Ar, P_el_Ar, _ = np.split(data_Ar, 3, axis = 1)
"""
#Load Pure ArCH4 93-7
data_ArCH4_937 = np.loadtxt('data/clusters_distributions/ArCH4-93-7-10bar_0.1x0.1x0.1mmCell_pi2.5GeV_normalized.txt' , skiprows = 1)
n_el_ArCH4_937, P_el_ArCH4_937, _ = np.split(data_ArCH4_937, 3, axis = 1)

"""
#Load experimental data for Ar
P_Ar_data = np.loadtxt('data/clusters_distributions/Ar_cluster_distribution_experimental.txt')/100

### DEGRAD ###
#Load DEGRAD data Pure Ar (from Carlos)
data_DEGRAD = np.loadtxt('data/clusters_distributions/Degrad_Pure_Ar_TEST.txt')
n_el_DEG, P_el_DEG = np.split(data_DEGRAD, 2, axis = 1)
P_el_DEG = P_el_DEG.flatten()
P_el_DEG_norm = P_el_DEG / sum(P_el_DEG)

#Load DEGRAD data Pure Ar v2 (from Carlos)
data_DEGRAD2 = np.loadtxt('data/clusters_distributions/Degrad_Pure_Ar_v3.txt')
n_el_DEG2, P_el_DEG2 = np.split(data_DEGRAD2, 2, axis = 1)
P_el_DEG2 = P_el_DEG2.flatten()
P_el_DEG_norm2 = P_el_DEG2 / sum(P_el_DEG2)
"""


#cut_for_extrapolation = len(P_Ar_data)

def PextrapolationAr(n):
    return 0.216 / n**2

def powerlaw(x,a):
    return a/x**2

def extrapolation2(data):
    opt, cov = curve_fit(powerlaw, n_el_DEG[cut_for_extrapolation-5:cut_for_extrapolation], data[-5:])
    return opt
    
"""
n = np.linspace(cut_for_extrapolation, e_cut, e_cut-cut_for_extrapolation)
P_Ar_extrapolated = PextrapolationAr(n)
P_Ar_final = np.concatenate((P_Ar_data, P_Ar_extrapolated))

n_el_DEG = np.linspace(cut_for_extrapolation, e_cut, e_cut-cut_for_extrapolation)
fit_params = extrapolation2(n)*4.8e-8
P_DEG_extrapolated = powerlaw(n_el_DEG, fit_params)
P_DEG_final = np.concatenate((P_el_DEG_norm[:cut_for_extrapolation], P_DEG_extrapolated))
"""

plt.figure()
plt.title('Ar/CF4 (99/1)', fontsize = 40)
#plt.plot(n_el_ArCH4_9010, P_el_ArCH4_9010, label = 'Ar/CH4 (90/10)')
#plt.plot(n_el_ArCF4_991, P_el_ArCF4_991, linestyle = 'dashdot' ,label = 'Ar/CF4 (99/1)')
plt.plot(n_el_ArCF4_991_le, P_el_ArCF4_991_le, label = '0.2 (GeV) $\\mu$')
plt.plot(n_el_ArCF4_991_he, P_el_ArCF4_991_he, label = '8.0 (GeV) $\\mu$')
plt.plot(n_el_ArCH4_937, P_el_ArCH4_937, label = 'ArCH4')

plt.xlabel('n electrons per cluster', fontsize = 25) ; plt.ylabel('P', fontsize = 25)
plt.xscale('log') ; plt.yscale('log')
plt.legend(fontsize = 18)
plt.xlim([0.9,1000])
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
"""
to_save = np.array((x, P_DEG_final))

np.savetxt('DEGRAD_ClusterSize.txt',to_save.T , delimiter='\t')
"""