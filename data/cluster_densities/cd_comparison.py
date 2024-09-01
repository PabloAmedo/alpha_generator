# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 10:09:28 2024

@author: diego
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

print('Running...')

def Momentum(energy_list, mass):
    p_list = []
    for energy in energy_list:
        p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)

P = 0.001 #bar 
fontsize_title = 25
fontsize_axis = 20
fontsize_legend = 15

#Load data ====================================================================
#HEED Cluster Density for Pions
PionMass = 139.57
PionCD = pd.read_csv('pre/ClusterDensity_N2_1.00__0.00_Pressure_7600.00Torr_pion_.txt', header=0, delimiter=' ')
PionCD['p/m (MeV/c)'] = PionCD['momentum'] * 1e-6 / PionMass

MuonMass = 105.66
MuonCD = pd.read_csv('pre/ClusterDensity_N2_1.00__0.00_Pressure_7600.00Torr_muon_.txt', header=0, delimiter=' ')
MuonCD['p/m (MeV/c)'] = MuonCD['momentum'] * 1e-6 / MuonMass

KaonMass = 493.7
KaonCD = pd.read_csv('pre/ClusterDensity_N2_1.00__0.00_Pressure_7600.00Torr_kaon_.txt', header=0, delimiter=' ')
KaonCD['p/m (MeV/c)'] = KaonCD['momentum'] * 1e-6 / KaonMass

ProtonMass = 938.27
ProtonCD = pd.read_csv('pre/ClusterDensity_N2_1.00__0.00_Pressure_7600.00Torr_proton_.txt', header=0, delimiter=' ')
ProtonCD['p/m (MeV/c)'] = ProtonCD['momentum'] * 1e-6 / ProtonMass

MuonCD_C3H8 = pd.read_csv('pre/ClusterDensity_C3H8_1.00__0.00_Pressure_7600.00Torr_muon_.txt', header=0, delimiter=' ')
MuonCD_C3H8['p/m (MeV/c)'] = MuonCD_C3H8['momentum'] * 1e-6 / MuonMass

#Merge
GeneralCD = pd.concat([PionCD, MuonCD, KaonCD, ProtonCD], ignore_index=True)

#Extrapolation
x = np.logspace(-2,3, 100000)
interpolator = interp1d(GeneralCD['p/m (MeV/c)'], GeneralCD['clusterDensity']/ 10 * P, kind = 'slinear', fill_value='extrapolate')
cd_extr_sli = interpolator(x)

interpolator_lin = interp1d(GeneralCD['p/m (MeV/c)'], GeneralCD['clusterDensity']/ 10 * P, kind = 'zero', fill_value='extrapolate')
interpolator_lin2 = interp1d(GeneralCD['p/m (MeV/c)'], GeneralCD['clusterDensity']/ 10 * P, kind = 'quadratic', fill_value='extrapolate')

cd_extr_lin = interpolator_lin(x)
cd_extr_lin2 = interpolator_lin2(x)
#From Santovetti et al.
#ItCD = pd.read_csv('../Ncl_Ar_good.csv', usecols=[0,1] ,names=['p/m (MeV/c)', 'clusterDensity'], delimiter=';')

plt.figure()
plt.plot(MuonCD['p/m (MeV/c)'], MuonCD['clusterDensity'] / 10 * P, 'bo', alpha = 0.5, label = '$\\mu$ @ N2')
plt.plot(ProtonCD['p/m (MeV/c)'], ProtonCD['clusterDensity'] / 10 * P, 'yo', alpha = 0.5, label = 'p @ N2')
plt.plot(MuonCD_C3H8['p/m (MeV/c)'], MuonCD_C3H8['clusterDensity'] / 10 * P, 'ro', alpha = 0.5, label = '$\\mu$ @ C3H8')
#plt.plot(PionCD['p/m (MeV/c)'], PionCD['clusterDensity'] / 10, 'ro', alpha = 0.5, label = 'HEED')
#plt.plot(KaonCD['p/m (MeV/c)'], KaonCD['clusterDensity'] / 10, 'go', alpha = 0.5, label = 'HEED')

#plt.plot(x, cd_extr_sli, 'k-', alpha = 0.5, label = 'slinear')
#plt.plot(x, cd_extr_lin, 'b-', alpha = 0.5, label = 'zero')
plt.plot(x, cd_extr_lin2, 'g-', alpha = 0.5, label = 'quadratic')
"""
plt.plot(500/105.66, 277.659, 'r<')
plt.plot(2500/105.66, 295.329, 'r<')
plt.plot(8000/105.66, 304.352, 'r<')
plt.plot(500/105.66, 256.4657803863417, 'b<')
plt.plot(2500/105.66, 289.91218103764504, 'b<')
plt.plot(8000/105.66, 302.69465649581804, 'b<')
"""
#plt.plot(ItCD['p/m (MeV/c)'], ItCD['clusterDensity'] * P, 'kx', label = 'Santovetti')
plt.xlabel('$\\beta \\gamma$', fontsize = fontsize_axis)
plt.ylabel('d$n_{ie}$/dx (1/cm)', fontsize = fontsize_axis)
#plt.xlim([0.5, 1000])
#plt.ylim([240, 500])
plt.grid()
plt.xlim([0.04, 75])
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize = fontsize_legend)