# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 14:11:07 2024

@author: diego
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from tools import *
from scipy.optimize import curve_fit

print('Running...')
#Create a function for this....................................................
os.chdir('../data_debug/')

#INPUTS =======================================================================
bkgcut1 = 0.35
bkgcut2 = 0.35
bkgcut3 = 0.35

base = 109
"""
THERE IS NO VERTICAL CUT IN THIS SCRIPT
"""
#Load simulated images ========================================================
sim50t1e1e = np.loadtxt('diff_2-7mm/sim50t_1e1e.txt')
sim50t1e10e = np.loadtxt('diff_2-7mm/sim50t_1e10e.txt')
sim50t1e100e = np.loadtxt('diff_2-7mm/sim50t_1e100e.txt')
sim50t1e1000e = np.loadtxt('diff_2-7mm/sim50t_1e1000e.txt')

sim100t1e1e = np.loadtxt('diff_2-7mm/sim100t_1e1e.txt')
sim100t1e10e = np.loadtxt('diff_2-7mm/sim100t_1e10e.txt')
sim100t1e100e = np.loadtxt('diff_2-7mm/sim100t_1e100e.txt')
sim100t1e1000e = np.loadtxt('diff_2-7mm/sim100t_1e1000e.txt')

sim500t1e10e = np.loadtxt('diff_2-7mm/sim500t_1e10e.txt')
sim500t1e100e = np.loadtxt('diff_2-7mm/sim500t_1e100e.txt')
sim500t1e1000e = np.loadtxt('diff_2-7mm/sim500t_1e1000e.txt')

sim1000t1e100e = np.loadtxt('diff_2-7mm/sim1000t_1e100e.txt')
sim1000t1e1000e = np.loadtxt('diff_2-7mm/sim1000t_1e1000e.txt')

sim10000t1e1000e = np.loadtxt('diff_2-7mm/sim10000t_1e1000e.txt')
#Load data ====================================================================
os.chdir('../tiffs/')
raw100ms = Image.open('July24/new set data/highgain_1structure/12x12/100ms/ss_single_1.tiff')
raw200ms = Image.open('July24/new set data/highgain_1structure/12x12/200ms/ss_single_1.tiff')
raw1000ms = Image.open('July24/new set data/highgain_1structure/12x12/700ms/ss_single_1.tiff')

img100ms_ = BaselineSubstraction(raw100ms, baseline = base)
img200ms_ = BaselineSubstraction(raw200ms, baseline = base)
img1000ms_ = BaselineSubstraction(raw1000ms, baseline = base)

img100ms = BkgSubstraction(img100ms_, bkg_cut = bkgcut1)
img200ms = BkgSubstraction(img200ms_, bkg_cut = bkgcut2)
img1000ms = BkgSubstraction(img1000ms_, bkg_cut = bkgcut3)

data100ms = np.array(img100ms)
data200ms = np.array(img200ms)
data1000ms = np.array(img1000ms)
# =============================================================================
#50 tracks
sum50t1e1e = np.sum(sim50t1e1e, axis = 0)
sum50t1e10e = np.sum(sim50t1e10e, axis = 0)
sum50t1e100e = np.sum(sim50t1e100e, axis = 0)
sum50t1e1000e = np.sum(sim50t1e1000e, axis = 0)
#100 tracks
sum100t1e1e = np.sum(sim100t1e1e, axis = 0)
sum100t1e10e = np.sum(sim100t1e10e, axis = 0)
sum100t1e100e = np.sum(sim100t1e100e, axis = 0)
sum100t1e1000e = np.sum(sim100t1e1000e, axis = 0)
#500 tracks
sum500t1e10e = np.sum(sim500t1e10e, axis = 0)
sum500t1e100e = np.sum(sim500t1e100e, axis = 0)
sum500t1e1000e = np.sum(sim500t1e1000e, axis = 0)
#1000 tracks
sum1000t1e100e = np.sum(sim1000t1e100e, axis = 0)
sum1000t1e1000e = np.sum(sim1000t1e1000e, axis = 0)
#10000 tracks
sum10000t1e1000e = np.sum(sim10000t1e1000e, axis = 0)

xsim = np.linspace(0,100,100)
"""
sum100ms = np.sum(data100ms, axis = 0)
sum200ms = np.sum(data200ms, axis = 0)
sum1000ms = np.sum(data1000ms, axis = 0)
"""
#I changed something so I rename data{}ms to don't have to change in advance
sum100ms = data100ms
sum200ms = data200ms
sum1000ms = data1000ms

cut100ms = sum100ms[np.argmax(sum100ms) - 50 : np.argmax(sum100ms) + 50]
cut200ms = sum200ms[np.argmax(sum200ms) - 50 : np.argmax(sum200ms) + 50]
cut1000ms = sum1000ms[np.argmax(sum1000ms) - 50 : np.argmax(sum1000ms) + 50]

xdata = np.linspace(0, len(cut100ms), len(cut100ms))
xdatareal = np.linspace(np.argmax(sum100ms) - 50, np.argmax(sum100ms) + 50, 100)

#Plots ========================================================================
#50 tracks --> 100 ms
plt.figure()
plt.title('50 tracks', fontsize = 25)
plt.plot(xsim, sum50t1e1e/max(sum50t1e1e), label = '1:1')
plt.plot(xsim - (np.argmax(sum50t1e10e) - np.argmax(sum50t1e1e)), sum50t1e10e/max(sum50t1e10e), label = '1:10')
plt.plot(xsim - (np.argmax(sum50t1e100e) - np.argmax(sum50t1e1e)), sum50t1e100e/max(sum50t1e100e), label = '1:100')
plt.plot(xsim - (np.argmax(sum50t1e1000e) - np.argmax(sum50t1e1e)), sum50t1e1000e/max(sum50t1e1000e), label = '1:1000')
plt.plot(xdata + (np.argmax(sum50t1e10e) - np.argmax(cut100ms)), cut100ms / max(cut100ms), label = 'Data', color = 'k', lw = 2)
plt.legend(fontsize = 18)
plt.xlabel('X (px)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize = 18)
#100 tracks --> 200 ms
plt.figure()
plt.title('100 tracks', fontsize = 25)
plt.plot(xsim, sum100t1e1e/max(sum100t1e1e), label = '1:1')
plt.plot(xsim - (np.argmax(sum100t1e10e) - np.argmax(sum100t1e1e)), sum100t1e10e/max(sum100t1e10e), label = '1:10')
plt.plot(xsim - (np.argmax(sum100t1e100e) - np.argmax(sum100t1e1e)) , sum100t1e100e/max(sum100t1e100e), label = '1:100')
plt.plot(xsim - (np.argmax(sum1000t1e100e) - np.argmax(sum100t1e1e)), sum100t1e1000e/max(sum100t1e1000e), label = '1:1000')
plt.plot(xdata + (np.argmax(sum50t1e10e) - np.argmax(cut200ms)), cut200ms / max(cut200ms), label = 'Data', color = 'k', lw = 2)
plt.legend(fontsize = 18)
plt.xlabel('X (px)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize = 18)
#500 tracks --> ~1000 ms
plt.figure()
plt.title('500 tracks', fontsize = 25)
plt.plot(xsim, sum500t1e10e/max(sum500t1e10e), label = '1:10')
plt.plot(xsim - (np.argmax(sum500t1e100e) - np.argmax(sum500t1e10e)), sum500t1e100e/max(sum500t1e100e), label = '1:100')
plt.plot(xsim - (np.argmax(sum500t1e1000e) - np.argmax(sum500t1e10e)), sum500t1e1000e/max(sum500t1e1000e), label = '1:1000')
plt.plot(xdata + (np.argmax(sum50t1e10e) - np.argmax(cut1000ms)), cut1000ms / max(cut1000ms), label = 'Data', color = 'k', lw = 2)
plt.legend(fontsize = 18)
plt.xlabel('X (px)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize = 18)

#Total phe comparison =========================================================
#50 tracks
plt.figure()
plt.title('50 tracks -- no normalized', fontsize = 25)
plt.plot(xsim, sum50t1e1e, label = '1:1')
plt.plot(xsim - (np.argmax(sum50t1e10e) - np.argmax(sum50t1e1e)), sum50t1e10e * 10, label = '1:10')
plt.plot(xsim - (np.argmax(sum50t1e100e) - np.argmax(sum50t1e1e)), sum50t1e100e * 100, label = '1:100')
plt.plot(xsim - (np.argmax(sum50t1e1000e) - np.argmax(sum50t1e1e)) +3 , sum50t1e1000e * 1000, label = '1:1000')
plt.legend(fontsize = 18)
plt.xlabel('X (px)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize = 18)
#100 tracks
plt.figure()
plt.title('100 tracks -- no normalized', fontsize = 25)
plt.plot(xsim, sum100t1e1e, label = '1:1')
plt.plot(xsim - (np.argmax(sum100t1e10e) - np.argmax(sum100t1e1e)), sum100t1e10e * 10, label = '1:10')
plt.plot(xsim - (np.argmax(sum100t1e100e) - np.argmax(sum100t1e1e)), sum100t1e100e * 100, label = '1:100')
plt.plot(xsim - (np.argmax(sum100t1e1000e) - np.argmax(sum100t1e1e)), sum100t1e1000e * 1000, label = '1:1000')
plt.legend(fontsize = 18)
plt.xlabel('X (px)', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize = 18)







