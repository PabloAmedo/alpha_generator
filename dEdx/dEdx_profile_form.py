# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:29 2024

@author: diego
"""
import os
os.chdir('../')
from Alpha_track_simulator import*

print('Running...\n')
#INPUTS========================================================================
gas          =       'ArCH4-90/10'
n_tracks     =       10000                                                     
Emuon        =       3000                                                      #(MeV) (CR)
mass         =       139.57                                                    #(MeV/c2)  pion
p            =       np.sqrt((Emuon + mass)**2 - mass**2)
dimensions   =       ([1.5,1.5,1.5])
rho_Ar       =       1.662                                                     #rho Ar/CH4 (93/7)
W_Ar         =       26.4e-3                                                   #keV
e_cut        =       50000 * 8                                                 #n_max electrons in a cluster
bins         =       30000 #FIX THIS

fact_data    =       1.0
peak_red_fact=       1.0

#==============================================================================
#==============================================================================

#Load Heed vs Data comparison
HEEDdata             =       np.loadtxt('data/HEED_v_data.csv', delimiter=';')
x_HEED, counts_HEED  =       np.split(HEEDdata, 2, axis=1 )
HEED_norm            =       counts_HEED / max(counts_HEED)


#Load Harris
Harrisdata               =       np.loadtxt('data/Harris72.csv', delimiter=';')
x_Harris, counts_Harris  =       np.split(Harrisdata, 2, axis=1 )
Harris_norm              =       counts_Harris / max(counts_Harris)
#LCP GENERATOR=================================================================
muons        =   []
muon         =   muon_generator(energy = Emuon, geometry = dimensions, gas = gas, mass = mass)          #First generate the muon generator object
muon.produce_muon(n = n_tracks, store = muons, line = True, e_cut = e_cut)

Sim_dEdx     =   []                                                             #Initialization list for simulated dE/dx -- truth
Det_dEdx     =   []                                                             #Initialization list for detected dE/dx  -- after avalanche effect

for track in muons:
    length = track.track_length
    e_per_cluster = track.n_e_cl

    Sim_dE = track.n_electrons * W_Ar
    Sim_dEdx.append(Sim_dE)
    
    Det_dE = sum(statistics_avalanche(e_per_cluster) * W_Ar)
    Det_dEdx.append(Det_dE / 1)


Sim_dEdx = np.array(Sim_dEdx)
Det_dEdx = np.array(Det_dEdx)

#Some calculations to compare the dEdx value: simulated vs known values
dEdx_simulated = np.mean(Sim_dEdx / rho_Ar)

hval, hbins = np.histogram(Det_dEdx   , bins = bins)

mpv_sim_aux = np.argmax(hval)
mpv_sim = hbins[mpv_sim_aux]

mpv_data = 2.2


Det_dEdx_norm = Det_dEdx / mpv_sim

e_per_cm = []
for muon in muons:
    e_per_cm.append(muon.n_electrons / muon.track_length)

e_cm = np.array(e_per_cm)
e_cm_prom = np.mean(e_cm)
print('<e/cm>:\t', e_cm_prom)
#PLOT==========================================================================

plt.figure()
plt.title('dE/dx profile (DETECTED)',   fontsize = 23)
hval, hbins,_ =plt.hist(Sim_dEdx ,    bins =  bins,    label = 'simulated',     color = 'r',    histtype = 'step')
#hval, hbins,_ =plt.hist(Sim_dEdx   ,   bins = bins,        label = 'simulated',    color = 'r',    histtype = 'step')
#plt.xscale('log')
plt.xlim([0 , 15])
plt.legend()


plt.figure()
#plt.plot(x_AC_data, ACdata, 'k.', label = 'Harris72-AC')
hval, hbins,_ =plt.hist(Det_dEdx ,    bins =  bins,    label = 'Sim',     color = 'k',    histtype = 'step')
plt.plot(x_HEED, (counts_HEED / max(counts_HEED)) * max(hval) * peak_red_fact, 'rv', label= 'HEED')
plt.plot(x_Harris, (counts_Harris / max(counts_Harris)) * max(hval)* peak_red_fact, 'bx', label= 'Harris')
plt.xlim([0 , 15])
plt.legend()

print("mpv_sim:\t", mpv_sim)
print("mpv_data:\t", mpv_data)
print("<dE/dx>:\t", np.mean(Det_dEdx) / 1.5 )