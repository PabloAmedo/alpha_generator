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
"""
#Load N2 10 bar 0.2 GeV
data_N2_02 = np.loadtxt('data/clusters_distributions/nitrogen-10bar_0.1x0.1x0.1mmCell_0.2GeVmuons_normalized.txt' , skiprows = 1)
n_el_N2_02, P_el_N2_02, _ = np.split(data_N2_02, 3, axis = 1)
"""
#Load N2 10 bar 2.5 GeV
data_N2_25 = np.loadtxt('data/clusters_distributions/nitrogen-10bar_0.1x0.1x0.1mmCell_2.5GeVmuons_normalized.txt' , skiprows = 1)
n_el_N2_25, P_el_N2_25, _ = np.split(data_N2_25, 3, axis = 1)

#Load N2 10 bar 8.0 GeV
data_N2_8 = np.loadtxt('data/clusters_distributions/propane-10bar_0.1x0.1x0.1mmCell_2.5GeVmuons_normalized.txt' , skiprows = 1)
n_el_N2_8, P_el_N2_8, _ = np.split(data_N2_8, 3, axis = 1)


# Crear la figura general
plt.figure(figsize=(12, 6))

# Primer subplot
plt.subplot(1, 2, 1)  # 1 fila, 2 columnas, primer subplot
plt.plot(n_el_N2_25, P_el_N2_25, label='N2')
plt.xlabel('n electrons per cluster', fontsize=25)
plt.ylabel('P', fontsize=25)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=18)
plt.xlim([0.9, 1000])
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()

# Segundo subplot
plt.subplot(1, 2, 2)  # 1 fila, 2 columnas, segundo subplot
plt.plot(n_el_N2_8, P_el_N2_8, label='C3H8')
plt.xlabel('n electrons per cluster', fontsize=25)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=18)
plt.xlim([0.9, 1000])
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()

# Ajustar el espacio entre subplots para que no se sobrepongan
plt.tight_layout()

# Mostrar el gráfico
plt.show()















"""
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

n_el_DEG = np.linspace(cut_for_extrapolation, e_cut, e_cut-cut_for_extrapolation)
fit_params = extrapolation2(n)*4.8e-8
P_DEG_extrapolated = powerlaw(n_el_DEG, fit_params)
P_DEG_final = np.concatenate((P_el_DEG_norm[:cut_for_extrapolation], P_DEG_extrapolated))


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

to_save = np.array((x, P_DEG_final))

np.savetxt('DEGRAD_ClusterSize.txt',to_save.T , delimiter='\t')
"""