# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:10:49 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

print('Running...')

def ClusterParametrizationAr(n):
    return 0.216/n**2 

def ClusterParametrizationCH4(n):
    return 0.119/n**2


e_cut = 10000

#Load Ar/CH4(90/10) data (from HEED presentation)
data_ArCH4_9010 = np.loadtxt('data/clusters_distributions/ArCH4_90-10_cluster_distribution_HEED.csv', delimiter=';')
_, P_el = np.split(data_ArCH4_9010, 2, axis = 1)
P_el = P_el.flatten()
n_el = np.linspace(1, 133, 133)
probabilities_133 = ClusterParametrizationAr(np.linspace(134, e_cut, e_cut-133))
probs_final = np.concatenate((P_el, probabilities_133))

#Load Ar/CF4(99/1) data (from HEED - Carlos)
data_ArCF4_9901 = np.loadtxt('data/clusters_distributions/ArCF4-99-1-10bar_0.1x0.1x0.1mmCell_2.5GeVmuons_normalized.txt' , skiprows = 1)
n_el, P_el, _ = np.split(data_ArCF4_9901, 3, axis = 1)

#Load experimental data for Ar
P_Ar_data = np.loadtxt('data/clusters_distributions/Ar_cluster_distribution_experimental.txt')/100

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



cut_for_extrapolation = len(P_Ar_data)

def PextrapolationAr(n):
    return 0.216 / n**2

def powerlaw(x,a):
    return a/x**2

def extrapolation2(data):
    opt, cov = curve_fit(powerlaw, n_el_DEG[cut_for_extrapolation-5:cut_for_extrapolation], data[-5:])
    return opt
    
n = np.linspace(cut_for_extrapolation, e_cut, e_cut-cut_for_extrapolation)
P_Ar_extrapolated = PextrapolationAr(n)
P_Ar_final = np.concatenate((P_Ar_data, P_Ar_extrapolated))
"""
n_el_DEG = np.linspace(cut_for_extrapolation, e_cut, e_cut-cut_for_extrapolation)
fit_params = extrapolation2(n)*4.8e-8
P_DEG_extrapolated = powerlaw(n_el_DEG, fit_params)
P_DEG_final = np.concatenate((P_el_DEG_norm[:cut_for_extrapolation], P_DEG_extrapolated))
"""
x = np.linspace(1,10000,10000)
x2 = np.linspace(1,1000,1000)

plt.figure()
plt.title('Cluster size distribution', fontsize = 20)
#plt.plot(x, probs_final, label = 'Ar/CH4 (90-10) presentation (HEED)')
plt.plot(x, P_el, label = 'Ar/CF4 (99-1) - HEED')
plt.plot(x, P_Ar_final, label = 'Fischle + extrapolation')
plt.plot(x[18], P_Ar_data[len(P_Ar_data)-1], 'ro')
plt.plot(x2, P_el_DEG_norm, 'kx', label = 'DEGRAD')
plt.plot(x2, P_el_DEG_norm2, 'bx', label = 'DEGRAD - 5000 eV cut')
plt.xlabel('n electrons per cluster', fontsize = 13) ; plt.ylabel('P', fontsize = 13)
plt.xscale('log') ; plt.yscale('log')
plt.legend()
plt.grid()
"""
to_save = np.array((x, P_DEG_final))

np.savetxt('DEGRAD_ClusterSize.txt',to_save.T , delimiter='\t')
"""