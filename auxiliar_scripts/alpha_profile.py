# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 09:37:38 2024

@author: diego
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import pandas as pd
from optical_gain_tools import *
import os

print('Running...')

#INPUTS =======================================================================
#Load EXP DATA
folder = '../tiffs/July24/new set data/highgain_1structure/12x12/700ms/'
path = folder + 'ss_single_1'
ext = '.tiff'
base = 109  #camera offset --> the number comes from previous analyses (see noise measurements)

#Load SIMULATED DATA
red = 10
h2 = np.loadtxt('../data_debug/diff_2-7mm/sim350t_1e{}e.txt'.format(red))
##-----------------------------------------------------------------------------
#Cut parameters
axis_cut = 'y'
other_axis = 'X' #This is only needed for the axis
y_offcenter = 10 #This has to be updated according to max (??)
x_offcenter = 8
xoff = -3

n    = 350    #Number of simulated tracks --> only needed for plot's title

#The next is all in pixel units, it has to be udpated to distance units!!
hpix = 62    #remove this num of pxs horizontally at each side --> useful to minimize tube reflections
vpix = 25    #height of the cut in pxs

#Display settings
show = True
save = False
cmap = 'jet'

#PRE-STEPS ====================================================================
"""
#This was used to iterate over a folder with a large amount of data
tiff_files = [f for f in os.listdir(folder) if f.endswith('.tiff')]
nfiles = 1
"""
gainlist = []


#CUSTOM CMAP --> allow to show 0s as white
cmap = plt.get_cmap(cmap)
norm = mcolors.Normalize(vmin=0, vmax=1)
colors = cmap(norm(np.linspace(0, 1, 256)))
new_colors = np.vstack((np.array([1, 1, 1, 1]), colors))

custom_cmap = mcolors.ListedColormap(new_colors)

#ITERATION ====================================================================
"""
#Uncomment in the case we want to run it for diff
for i in range(1, nfiles):
   image_path100 = path + str(i) + ext
"""
image_path100 = path + ext

image100 = Image.open(image_path100)
width100, height100 = image100.size

image100 = BaselineSubstraction(image100, baseline = base)
imageArray100, aux100 = ImageCut([0, height100/2, width100, height100], 
                                 image100, vpix = vpix, hpix= hpix, 
                                 show = show, axis_cut = axis_cut, offcenter = [x_offcenter, y_offcenter])
imageArray100 = BkgSubstraction(imageArray100, 0.35)
data = LoadOpticalParams(12)



x100 = np.linspace(0 + xoff, len(imageArray100) +xoff, len(imageArray100))
#==============================================================================

h2cut = h2[:][26:76]
h2sum = np.sum(h2cut, axis = 0)

nedata = sum(imageArray100) / float(data['qeff']) / float(data['geomeff']) / float(data['T'])
nesim = sum(h2sum) * red
gain = nedata / nesim
gainlist.append(gain)

if show == True:
    """
    #PLOT 1: Data spectrum 
    plt.figure()
    plt.title('Data spectrum', fontsize = 25)
    plt.plot(x100, imageArray100, 'b', label = '100 ms')
    plt.legend(fontsize = 18)
    plt.xlabel('{} (px)'.format(other_axis), fontsize = 20)
    plt.ylabel('Counts', fontsize = 20)
    plt.tick_params(axis='both', which='major', labelsize = 18)
    
    #PLOT 2: Sim spectrum 
    plt.figure()
    plt.title('Simulated scpectrum', fontsize = 25)
    plt.plot(x100, h2sum, 'b')
    plt.legend(fontsize = 18)
    plt.xlabel('{} (px)'.format(other_axis), fontsize = 20)
    plt.ylabel('Counts', fontsize = 20)
    plt.tick_params(axis='both', which='major', labelsize = 18)
    """
    #PLOT 3: Comparison spectrum between data and sim
    plt.figure()
    plt.title('{:.2f} tracks'.format(n), fontsize = 25)
    plt.plot(np.linspace(0, 100, 100), h2sum/max(h2sum), label = 'sim')
    plt.plot(x100, imageArray100/max(imageArray100), 'k', label = 'data')
    plt.tick_params(axis='both', which='major', labelsize = 18)
    plt.xlabel('X (px)', fontsize = 20)
    plt.ylabel('Counts', fontsize = 20)
    plt.legend(fontsize = 18)
    #RESULTS: Numerically results 
    print('\nNe data:\t{:.3e}'.format(nedata))
    print('Ne sim:\t\t{:.3e}'.format(nesim))
    print('\nGain:\t\t{:.3f}'.format(gain))
    
print('Average Optical Gain:\t{:.3f}'.format(np.mean(gainlist)))
        