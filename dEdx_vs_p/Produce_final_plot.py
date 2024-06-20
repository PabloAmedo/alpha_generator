# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:36:18 2024

@author: usuario
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.optimize import curve_fit
import warnings

import os

warnings.simplefilter(action="ignore")

plt.close('all')

print('Running...')

#Inputs
binsx = int(150*2)
binsy = int(75*2)
lower_p_cut = 0

fontsize_title = 25
fontsize_subtitle = 20
fontsize_axis = 20
fontsize_subaxis = 15
fontsize_legend = 15

pcut1_lower = 70
pcut1_upper = 100

bins1       = 35
bins2       = 75

#Load Data into df
dataset = pd.read_csv('simulated_data/PEP4_TEST.txt',     delimiter='\t', header=10)
dataset_ = pd.read_csv('simulated_data/PEP4_more_stat.txt',     delimiter='\t', header=10)
dataset1 = pd.read_csv('simulated_data/ArCF4_99-1_PLOT_2_dEdx-v-p.txt',     delimiter='\t', header=10)
dataset2 = pd.read_csv('simulated_data/ArCF4_99-1_PLOT_3_dEdx-v-p.txt',     delimiter='\t', header=10)
dataset3 = pd.read_csv('simulated_data/ArCF4_99-1_PLOT_dEdx-v-p.txt',       delimiter='\t', header=10)
dataset4 = pd.read_csv('simulated_data/ArCF4_99-1_PLOT_2_dEdx-v-p_lt.txt',  delimiter='\t', header=10)
dataset5 = pd.read_csv('simulated_data/ArCF4_99-1_PLOT_dEdx-v-p_lt.txt',    delimiter='\t', header=10)
dataset6 = pd.read_csv('simulated_data/ArCF4_99-1_PLOT_3_dEdx-v-p_lt.txt',  delimiter='\t', header=10)
#Merge into the same dataframe
data_ = pd.concat((dataset1, dataset2, dataset3, dataset4, dataset5, dataset6), axis = 0, ignore_index=True)
data_ = pd.concat((dataset, dataset_), axis = 0, ignore_index= True)
data = data_[(data_['P'] >= lower_p_cut)]

mass_ID = {
    0.511:4,     #electron
    105.66:0,   #muon
    139.57:1,   #pion
    493.7:2,    #kaon
    938.27:3    #proton
    }

data['ID'] = data['Mass'].map(mass_ID)

data_muon_ = data[(data['ID'] == 0)].reset_index(drop=True)
data_pion_ = data[(data['ID'] == 1)].reset_index(drop=True)
data_kaon_ = data[(data['ID'] == 2)].reset_index(drop=True)
data_proton_ = data[(data['ID'] == 3)].reset_index(drop=True)
data_electron_ = data[(data['ID'] == 4)].reset_index(drop=True)

#Load momentum list and momentum resolution
p_list_muon = np.loadtxt('../data/p_lists/p_muon.txt').round(5)
p_list_pion = np.loadtxt('../data/p_lists/p_pion.txt').round(5)
p_list_kaon = np.loadtxt('../data/p_lists/p_kaon.txt').round(5)
p_list_proton = np.loadtxt('../data/p_lists/p_proton.txt').round(5)
p_list_electron = np.loadtxt('../data/p_lists/p_electron.txt').round(5)

p_res_muon = np.loadtxt('../data/p_resolution/p_muon_res')/100
p_res_pion = np.loadtxt('../data/p_resolution/p_pion_res')/100
p_res_kaon = np.loadtxt('../data/p_resolution/p_kaon_res')/100
p_res_proton = np.loadtxt('../data/p_resolution/p_proton_res')/100
p_res_electron = np.loadtxt('../data/p_resolution/p_electron_res')/100

p_muon_ = np.array([p_list_muon, p_res_muon]).T
p_pion_ = np.array([p_list_pion, p_res_pion]).T
p_kaon_ = np.array([p_list_kaon, p_res_kaon]).T
p_proton_ = np.array([p_list_proton, p_res_proton]).T
p_electron_ = np.array([p_list_electron, p_res_electron]).T

p_muon = pd.DataFrame(p_muon_, columns=['P_true','sigmaP'])
p_pion = pd.DataFrame(p_pion_, columns=['P_true','sigmaP'])
p_kaon = pd.DataFrame(p_kaon_, columns=['P_true','sigmaP'])
p_proton = pd.DataFrame(p_proton_, columns=['P_true','sigmaP'])
p_electron = pd.DataFrame(p_electron_, columns=['P_true','sigmaP'])

#Merge -----> dropna() to avoid the missmatch between the decimals in p_true
data_muon = pd.merge(data_muon_,p_muon, on = 'P_true', how='left').dropna()
data_pion = pd.merge(data_pion_, p_pion, on = 'P_true', how='left').dropna()
data_kaon = pd.merge(data_kaon_, p_kaon, on = 'P_true', how='left').dropna()
data_proton = pd.merge(data_proton_, p_proton, on = 'P_true', how='left').dropna()
data_electron = pd.merge(data_electron_, p_electron, on = 'P_true', how='left').dropna()
#Mergeing data from all particles
data_final = pd.concat((data_muon, data_pion, data_kaon, data_proton, data_electron))

#Applying the resolution in momentum !!!!!!!!! FIXME --> this can be done once and then load all data
data_final['P_good'] = np.random.normal(loc = data_final['P_true'], 
                                        scale = np.sqrt((data_final['P_true']**2 * 0.035)))

Cut1 = data_final[(data_final['P_good'] > pcut1_lower) & (data_final['P_good'] < pcut1_upper)]
#PLOT==========================================================================

x_space = np.logspace(1.75,5.05, binsx)                                            #FIXME -- this is on MeV/c
y_space = np.linspace(min(data['dEdx_trunc']) - 0.05 * min(data['dEdx_trunc']), max(data['dEdx_trunc']) + 0.05 * max(data['dEdx_trunc']), binsy)

plt.figure()
plt.title('{:.3e} events'.format(len(data_final)), fontsize = 20)
plt.hist2d(data_final['P_good'], data_final['dEdx_trunc'],bins=(x_space, y_space), cmin=1,cmap = "jet")
#plt.hist2d(data_final['P_true'], data_final['dEdx_trunc'],bins=(x_space, y_space), cmin=1,cmap = "Oranges"  , alpha = 0.5)
plt.xlabel('P (MeV/c)', fontsize = 20 )
plt.ylabel('dE/dx (keV/cm)', fontsize = 20 )
plt.axvspan(pcut1_lower, pcut1_upper, color='grey', alpha=0.5)
plt.xscale('log')
plt.colorbar()

#DOING THE CUT
#Cut1
Cut1_muon = Cut1[Cut1['ID'] == 0]
Cut1_pion = Cut1[Cut1['ID'] == 1]
Cut1_kaon = Cut1[Cut1['ID'] == 2]
Cut1_proton = Cut1[Cut1['ID'] == 3]
Cut1_electron = Cut1[Cut1['ID'] == 4]


fig, axs = plt.subplots(1,2)

axs[0].plot(Cut1['P_good'], Cut1['dEdx_trunc'], 'ko', alpha = 0.3)
axs[0].set_xlabel('P (MeV/c)', fontsize = fontsize_subaxis)
axs[0].set_ylabel('dE/dx (keV/cm)', fontsize = fontsize_subaxis)
axs[0].tick_params(axis='both', which='major', labelsize=fontsize_legend)

#axs[1].hist(Cut1['dEdx_trunc'], bins = 150, color='k', histtype = 'step', label = 'total', zorder = 3, linewidth = 1.5)
MuonCutData = axs[1].hist(Cut1_muon['dEdx_trunc'], bins = bins1, histtype = 'step', label = 'muon')
PionCutData = axs[1].hist(Cut1_pion['dEdx_trunc'], bins = bins1, histtype = 'step', label = 'pion')
#KaonCutData = axs[1].hist(Cut1_kaon['dEdx_trunc'], bins = bins1, histtype = 'step', label = 'kaon')
#ProtonCutData = axs[1].hist(Cut1_proton['dEdx_trunc'], bins = bins1, histtype = 'step', label = 'proton')
ElectronCutData = axs[1].hist(Cut1_electron['dEdx_trunc'], bins = bins1, histtype = 'step', label = 'electron')
#axs[1].hist(Cut1_kaon['dEdx_trunc'], bins = 75, histtype = 'step', label = 'kaon')
axs[1].set_xlabel('dE/dx (keV/cm)', fontsize = fontsize_subaxis)
axs[1].set_ylabel('Counts', fontsize = fontsize_subaxis)

fig.tight_layout()

# dEdx RESOLUTION =============================================================
def gaussian(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

xToMuon = np.linspace(min(MuonCutData[1]) - 0.25, max(MuonCutData[1]) + 0.25, 1000)
optMuon, errMuon = curve_fit(gaussian, MuonCutData[1][:-1], MuonCutData[0], p0=[max(MuonCutData[0]), np.mean(MuonCutData[1]), 1])
axs[1].plot(xToMuon, gaussian(xToMuon, *optMuon), 'b-', label = 'dEdx: {:.2f} (%)'.format(abs(optMuon[2]/optMuon[1]*100)), linewidth = 2, zorder = 5)

xToPion = np.linspace(min(PionCutData[1]) - 0.25, max(PionCutData[1]) + 0.25, 1000)
optPion, errPion = curve_fit(gaussian, PionCutData[1][:-1], PionCutData[0], p0=[max(PionCutData[0]), np.mean(PionCutData[1]), 1])
axs[1].plot(xToPion, gaussian(xToPion, *optPion), 'r-', label = 'dEdx: {:.2f} (%)'.format(abs(optPion[2]/optPion[1]*100)), linewidth = 2, zorder = 5)

xToElectron = np.linspace(min(ElectronCutData[1]) - 0.25, max(ElectronCutData[1]) + 0.25, 1000)
optElectron, errElectron = curve_fit(gaussian, ElectronCutData[1][:-1], ElectronCutData[0], p0=[max(ElectronCutData[0]), np.mean(ElectronCutData[1]), 1])
axs[1].plot(xToElectron, gaussian(xToElectron, *optElectron), 'g-', label = 'dEdx: {:.2f} (%)'.format(abs(optElectron[2]/optElectron[1]*100)), linewidth = 2, zorder = 5)

"""
xToKaon = np.linspace(min(KaonCutData[1]) - 0.25, max(KaonCutData[1]) + 0.25, 1000)
optKaon, errKaon = curve_fit(gaussian, KaonCutData[1][:-1], KaonCutData[0], p0=[max(KaonCutData[0]), np.mean(KaonCutData[1]), 1])
axs[1].plot(xToKaon, gaussian(xToKaon, *optKaon), label = 'dEdx: {:.2f} (%)'.format(abs(optKaon[2]/optKaon[1]*100)), linewidth = 2, zorder = 5)

xToProton = np.linspace(min(ProtonCutData[1]) - 0.25, max(ProtonCutData[1]) + 0.25, 1000)
optProton, errProton = curve_fit(gaussian, ProtonCutData[1][:-1], ProtonCutData[0], p0=[max(ProtonCutData[0]), np.mean(ProtonCutData[1]), 1])
axs[1].plot(xToProton, gaussian(xToProton, *optProton), label = 'dEdx: {:.2f} (%)'.format(abs(optProton[2]/optProton[1]*100)), linewidth = 2, zorder = 5)
"""
axs[1].legend(ncol =2, fontsize = fontsize_legend)
axs[1].tick_params(axis='both', which='major', labelsize=fontsize_legend)
#axs[1].set_xlim([min(np.concatenate((MuonCutData[1], PionCutData[1])))- 0.5,max(np.concatenate((MuonCutData[1], PionCutData[1]))) + 0.5])
fig.show()

#print('dEdx resolution\n')
#print('dEdx res mu:\t{:.4f} %'.format(abs(optMuon[2]/optMuon[1]*100)))
#print('dEdx res pi:\t{:.4f} %'.format(abs(optPion[2]/optPion[1]*100)))
#print('dEdx res k:\t\t{:.4f} %'.format(abs(optKaon[2]/optKaon[1]*100)))
#print('dEdx res p:\t\t{:.4f} %'.format(abs(optProton[2]/optProton[1]*100)))
