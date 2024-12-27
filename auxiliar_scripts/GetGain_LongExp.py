# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:34:44 2024

@author: diego
"""

import numpy as np
import os
import time
from PIL import Image
from PIL.ExifTags import TAGS

os.chdir('../')
from Alpha_track_simulator import *

st = time.time()


def GetMetadata(file, key = 'ImageDescription', path = None):
    with Image.open(file) as img:
        meta_dict = {TAGS[key]: img.tag[key] for key in img.tag_v2 if key in TAGS}
    if key == 'ImageDescription':
        meta = meta_dict[key][0]
        print(meta)
        return meta
    else:
        return meta_dict[key]

def OQ2Load(file, path = None):
    """
    Loads a file from OQ2 native format and returns a numpy array.
    
    OQ2 saves the images in .tif format as int16 which numpy is unable to read 
    so it has to be loaded using PIL.Image and then convert it to numpy.array
    """
    
    fullpath = path + file if path != None else file
    img = Image.open(fullpath)
    img_array = np.array(img).astype(np.int64) #maybe float here
    return img_array

def GetGain(path, primaries_profile = None ,offset = 0, QE = 0.7, GE = 1.6511e-4 , T = 1, inf = 400, sup = 650):
    gains = [] #gain using the max
    gains2 = [] #gain using the integral
    profs = []
    os.chdir(path)
    files = [file for file in os.listdir() if os.path.isfile(file)]
    
    for file in files:
        print(file)
        data = OQ2Load(file)
        data_prof = np.sum(data, axis = 0) / QE / GE / T - offset
        profs.append(data_prof)
        gains.append(max(data_prof) / max(primaries_profile))
        gains2.append(np.sum(data_prof[inf:sup]) / np.sum(primaries_profile[inf:sup]))
    return data_prof, gains, gains2, np.array(profs)


#INPUTS =======================================================================
n_tracks    =   875
red         =   100
native_bin  =   4
P           =   2.5

diff        =   0.22 

bins        =   [int(4096 / native_bin),int(2304 / native_bin)]                #Change to binning --> this is the size for a 12x12 reb
pixel_size  =   (4.6e-4 * native_bin * 19.48 )#* 0.85)                                  #real size * rebinning * magnification 
QE          =   0.7
GE          =   1.6511e-4
T           =   1

path        =   '../../NOV24 - EXPERIMENTAL CAMPAIGN/19122024/5bar/Both/500ms/'

x_range     =   [-(pixel_size * bins[0])/2, (pixel_size * bins[0])/2]
y_range     =   [-(pixel_size * bins[1])/2, (pixel_size * bins[1])/2]

save        =   0

cut1 = 228
cut2 = 727

print(os.getcwd())
#SIMULATE EXPERIMENTAL CONDITIONS  ============================================
source = Source(radius = 0.35, red_fact = red, rate = 1750, P = P)
track_list = []

exposition_time = n_tracks / (source.rate)
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg")

diff_handler=Diffusion_handler(sigma_diff = diff, sigma_PSF = 0)
diff_handler.diffuse(track_list)

for i in track_list: i.fill(diff = True)

noise = Noise(0.06)
#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins": bins, "range":[x_range, y_range]}, 
                   QE = QE, GE = GE, Tp = T, 
                   gain = 1)

hist, ecut = image2d.plot_x(variable_cut = 'x')#, min_var = -0.648, max_var = 0.378)
primaries = image2d.Hist2D_e
primaries_profile = np.sum(primaries, axis = 0) * red

recorded = image2d.Hist2D * red

#LOAD DATA ---------> 20-11
data_prof, gains, gains2, profs = GetGain(path, 
                           primaries_profile = primaries_profile,
                           offset = 102095958 * 0.99,#4720471,
                           QE = QE, GE = GE, T = T,
                           inf = cut1, sup = cut2)

# Plot distribution of the gains for each image
plt.figure()
#plt.title('28/11/2021 -- 2bar', fontsize = 15)
plt.hist(gains, bins = 25, histtype = 'step', label = 'Peaks')
plt.hist(gains2, bins = 25, histtype = 'step', label = 'Total phe')
plt.xlabel('Optical Gain (phe/e)', fontsize = 15)
plt.legend()
#plt.savefig('../../NOV24 - EXPERIMENTAL CAMPAIGN/OpticalGain-20112024.png')

x = np.linspace(0, len(data_prof), len(data_prof))

sec = image2d.Hist2D * red# / QE / GE / T 
sec_profile = np.sum(sec, axis = 0)
meanprofs = np.mean(profs, axis = 0)

# Plot profile comp
"""
plt.figure()
plt.plot(x, profs / max(profs),  label = 'Data')
plt.plot(x - 7, primaries_profile / max(primaries_profile), label = 'Sim')
plt.xlabel('px', fontsize = 15)
plt.ylabel('phe', fontsize = 15)
plt.xlim([cut1, cut2])
plt.legend()
"""

plt.figure()
plt.plot(x, meanprofs ,  label = 'Data')
plt.plot(x - 7, primaries_profile * QE* GE*T, label = 'Sim')
plt.xlabel('px', fontsize = 15)
plt.ylabel('phe', fontsize = 15)
plt.xlim([cut1, cut2])
plt.legend()




plt.figure()
plt.plot(x, meanprofs, label = 'Data')

plt.figure()
plt.plot(meanprofs/max(meanprofs))
plt.plot(x - 7, primaries_profile / max(primaries_profile))


print('='*50)
print('Averaged over:\t\t\t{:.3f} files'.format(len(gains)))
print('Aerage max:\t\t\t\t{:.3f} (phe/e)'.format(np.mean(gains)))
print('Aerage int:\t\t\t\t{:.3f} (phe/e)'.format(np.mean(gains2)))

print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(max(gains)))
print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(min(gains)))


"""
#LOAD DATA ---------> 21-11
data_prof, gains = GetGain('../../NOV24 - EXPERIMENTAL CAMPAIGN/21112024/200ms/', 
                           primaries_profile = primaries_profile,
                           offset = 0,
                           QE = 0.7, GE = 1.6511e-4, T = 1)
    

# Plot distribution of the gains for each image
plt.figure()
plt.title('21/11/2021 -- Max Gain', fontsize = 15)
plt.hist(gains, bins = 25, histtype = 'step')
plt.xlabel('Optical Gain (phe/e)', fontsize = 15)
plt.savefig('../../NOV24 - EXPERIMENTAL CAMPAIGN/Plots/OpticalGain-21112024.png')


print('='*15, ' 21/11/2024 ', '='*15)
print('Averaged over:\t\t\t{:.3f} files'.format(len(gains)))
print('Aerage:\t\t\t\t\t{:.3f} (phe/e)'.format(np.mean(gains)))
print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(max(gains)))
print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(min(gains)))

#LOAD DATA  ---------> 22-11
data_prof, gains = GetGain('../../NOV24 - EXPERIMENTAL CAMPAIGN/22112024/Acrylic only/200ms/4x4/', 
                           primaries_profile = primaries_profile,
                           offset = 0,
                           QE = 0.7, GE = 1.6511e-4, T = 0.5)
    

# Plot distribution of the gains for each image
plt.figure()
plt.title('22/11/2021 -- Acrylic only', fontsize = 15)
plt.hist(gains, bins = 25, histtype = 'step')
plt.xlabel('Optical Gain (phe/e)', fontsize = 15)
plt.savefig('../../NOV24 - EXPERIMENTAL CAMPAIGN/Plots/OpticalGain-22112024.png')
    
print('='*15, ' 22/11/2024 ', '='*15)
print('Averaged over:\t\t\t{:.3f} files'.format(len(gains)))
print('Aerage:\t\t\t\t\t{:.3f} (phe/e)'.format(np.mean(gains)))
print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(max(gains)))
print('Max Optical Gain:\t\t{:.3f} (phe/e)'.format(min(gains)))
"""
end = time.time()

print('\nElapsed time:\t{} s'.format(end - st))


