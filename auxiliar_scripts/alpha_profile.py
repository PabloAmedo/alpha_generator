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
from general_tools import *
import os

print('Running...')
plt.close('all')


#INPUTS =======================================================================
#Load SIM and DATA

ms = 350
red = 10

px_to_cm = 9.467

#DATA
folder = '../../Gain anlysis/Data/20240719/acryliconly/c-1136_r-128_acground_aa1910_fc2210_fa2210/100ms/'#+ str(int(ms/3)) + 'ms/'
#folder = '../tiffs/test/bin 12x12 100ms/'
path = folder + 'ss_single_'
ext = '.tiff'

#Create a file to save the data
save = open(folder + 'OpticalGain.txt', 'w')
save.write('INPUTS\n')
save.write('Exposure time:\t\t{}\n'.format(str(int(ms/3.5))))
save.write('Thinning factor:\t{}\n\n'.format(red))

base = 414  #camera offset --> the number comes from previous analyses (see noise measurements)

#SIM
h2 = np.loadtxt('../data_debug/diff_fixed_2-5mm/sim' + str(int(ms/(2))) + 't_1e{}e.txt'.format(red))
#h2 = np.loadtxt('../data_debug/diff_2-7mm/sim50t_1e{}e.txt'.format(red))
##-----------------------------------------------------------------------------
#Cut parameters
axis_cut = 'y'
other_axis = 'X' #This is only needed for the axis
y_offcenter = 15 #This has to be updated according to max (??)
x_offcenter = 5
xoff = 0
#bkgcut = 0.5

n    = 175    #Number of simulated tracks --> only needed for plot's title

#The next is all in pixel units, it has to be udpated to distance units!!
hpix = 55    #remove this num of pxs horizontally at each side --> useful to minimize tube reflections
vpix = 40    #25 #height of the cut in pxs

#Display settings
show = False
#save = False
cmap = 'RdBu'

#PRE-STEPS ====================================================================

#This was used to iterate over a folder with more than 1 img
tiff_files = [f for f in os.listdir(folder) if f.endswith('.tiff')]
nfiles = len(tiff_files)

gainlist = []           #Initialize a list where the gain values will be stored

#CUSTOM CMAP --> allow to show 0s as white -----------> Move to general_tools
cmap = plt.get_cmap(cmap)
norm = mcolors.Normalize(vmin=0, vmax=1)
colors = cmap(norm(np.linspace(0, 1, 256)))
new_colors = np.vstack((np.array([1, 1, 1, 1]), colors))

custom_cmap = mcolors.ListedColormap(new_colors)

#FLUCTUATIONS

#e- fluctuations
rows = len(h2)
h2ef_ = np.zeros_like(h2)
for i in range(rows):
    for j in range(len(h2[i])):
        acum = 0
        for _ in range(int(h2[i][j])):
            acum = acum + np.random.exponential(scale = 1 )
        h2ef_[i][j] = acum
        """
        else:
            h2ef_[i][j] = 0
        """
        

#ph fluctuations
h2ef = np.random.poisson(h2ef_)



#ITERATION ====================================================================

#nfiles = 2
#Uncomment in the case we want to run it for diff
for i in range(1, nfiles+1):
   data_path = path + str(i) + ext
   #image_path100 = path + ext
    
   data_image = Image.open(data_path)
   
   width100, height100 = data_image.size
    
   data_image_ = BaselineSubstraction(data_image, baseline = base)
   """
   data_profile, aux100 = ImageCut([0, height100/2, width100, height100], 
                                     data_image, vpix = vpix, hpix= hpix, 
                                     show = show, axis_cut = axis_cut, offcenter = [x_offcenter, y_offcenter], cmap = cmap)
   """

   aux_ = data_image_[44:124][:]
   aux = np.array([x[60:154] for x in aux_])
   data_profile = np.sum(aux, axis=0)
   #data_profile = BkgSubstraction(data_profile, bkgcut, show=show)
   data_image = np.array(data_image)
   opt_par = LoadOpticalParams(12)
   
   data_profile = data_profile - min(data_profile)
   
   #==============================================================================
   #h2cut = h2[:][26:76]
   h2cut_ = h2ef[44:124][:]
   h2cut = np.array([x[60:174] for x in h2cut_])
   h2sum = np.sum(h2cut, axis = 0)
    
   x = np.linspace(0 + xoff, len(data_profile) +xoff, len(data_profile))
   
   nedata = sum(data_profile) * 1.32 / float(opt_par['qeff'].iloc[0]) / float(opt_par['geomeff'].iloc[0]) / float(opt_par['T'].iloc[0])
   nesim = sum(h2sum) * red
   gain = nedata / nesim        ;        save.write(data_path[-16:] + '\t{}\n'.format(gain))
   
   gainlist.append(gain)
    
   if show == True:
    
       plt.figure(figsize = (10,7))
       plt.title('Data', fontsize = 25)
       a = plt.imshow(data_image_, cmap = cmap)
       plt.colorbar(a)
       plt.vlines(60, 44, 124, colors='r')
       plt.vlines(174, 44, 124, colors='r')
       plt.hlines(44, 60, 174, colors='r')
       plt.hlines(124, 60, 174, colors='r')
       plt.tick_params(axis='both', which='major', labelsize = 18)
       plt.xlabel('X (px)', fontsize = 20)
       plt.ylabel('Y (px)', fontsize = 20)
       
       
       plt.figure(figsize = (10,7))
       plt.title('Simulated', fontsize = 25)
       b = plt.imshow(h2ef, cmap=cmap)
       plt.colorbar(b)
       plt.vlines(60, 44, 124, colors='r')
       plt.vlines(174, 44, 124, colors='r')
       plt.hlines(44, 60, 174, colors='r')
       plt.hlines(124, 60, 174, colors='r')
       plt.tick_params(axis='both', which='major', labelsize = 18)
       plt.xlabel('X (px)', fontsize = 20)
       plt.ylabel('Y (px)', fontsize = 20)
       
       #PLOT 1: Data-cut spectrum 
       plt.figure(figsize = (12,7))
       plt.title('Data spectrum', fontsize = 25)
       plt.plot(x, data_profile, 'b', label = '100 ms')
       plt.legend(fontsize = 18)
       plt.xlabel('{} (px)'.format(other_axis), fontsize = 20)
       plt.ylabel('Counts', fontsize = 20)
       plt.tick_params(axis='both', which='major', labelsize = 18)
        
       #PLOT 2: Sim spectrum 
       plt.figure(figsize = (12,7))
       plt.title('Simulated scpectrum', fontsize = 25)
       plt.plot(np.linspace(0 + xoff, len(h2sum) +xoff, len(h2sum)), h2sum, 'b')
       plt.legend(fontsize = 18)
       plt.xlabel('{} (px)'.format(other_axis), fontsize = 20)
       plt.ylabel('Counts', fontsize = 20)
       plt.tick_params(axis='both', which='major', labelsize = 18)
       
       #PLOT 3: Comparison spectrum between data and sim
       plt.figure(figsize = (12,7))
       plt.title('{:.2f} tracks'.format(n), fontsize = 25)
       plt.plot(np.linspace(0, 114, 114), h2sum/max(h2sum), label = 'sim')
       plt.plot(x -5, data_profile / max(data_profile), 'k', label = 'data')
       plt.tick_params(axis='both', which='major', labelsize = 18)
       plt.xlabel('X (px)', fontsize = 20)
       plt.ylabel('Counts', fontsize = 20)
       plt.legend(fontsize = 18)
       
       
       #RESULTS: Numerically results 
       print('\nNe data:\t{:.3e}'.format(nedata))
       print('Ne sim:\t\t{:.3e}'.format(nesim))
       print('Gain:\t\t{:.3f}'.format(gain))
       print('-'*25)
    
print('Average Optical Gain:\t{:.3f}'.format(np.mean(gainlist)))
save.write('\nAverage gain:\t\t{}'.format(np.mean(gainlist)))
save.close()


x = np.linspace(0, len(h2[0]), len(h2[0]))
plt.figure()
plt.plot(x, np.sum(h2, axis = 0), label= 'h2')
plt.plot(x, np.sum(h2ef_, axis = 0), label= 'h2ef_')
plt.plot(x, np.sum(h2ef, axis = 0), label= 'h2ef')