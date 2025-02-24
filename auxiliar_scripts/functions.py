# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:17:46 2025

@author: diego
"""

import numpy as np
import os
import time
from PIL import Image
from PIL.ExifTags import TAGS
import matplotlib.pyplot as plt

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

def GetGain(path, primaries_profile, cuts, QE = 0.7, GE = 1.6511e-4 , T = 1, show = False):
    gains_peaks = [] #gain using the max
    gains_integral = [] #gain using the integral
    profs = []
    y1, y2 = cuts[0]
    x1, x2 = cuts[1]
    
    os.chdir(path)
    files = [file for file in os.listdir() if os.path.isfile(file)]
    
    for file in files:
        #print(file)
        data = OQ2Load(file)
        data_prof = np.sum(data[y1:y2], axis = 0) / QE / GE / T
        base = np.mean(data_prof[355:365])
        data_prof = data_prof - base
        profs.append(data_prof)
        gains_peaks.append(max(data_prof) / max(primaries_profile))
        gains_integral.append(np.sum(data_prof[x1:x2]) / np.sum(primaries_profile[x1:x2]))
    if show:
        plt.figure()
        plt.title('Cut')
        plt.imshow(data[y1:y2])#[:,300:700])
        plt.colorbar()
        
        plt.figure()
        plt.title('Raw')
        plt.imshow(data)
        plt.colorbar()
    
    return data_prof, gains_peaks, gains_integral, np.array(profs)