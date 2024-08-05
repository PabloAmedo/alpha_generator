# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:18:10 2024

@author: diego
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from scipy.optimize import curve_fit
import matplotlib.patches as patches
import matplotlib.colors as mcolors


def BaselineSubstraction(image, baseline = 100):                               #This is good, maybe no need a func
    image = np.array(image)
    image[image < baseline] = 0
    image[image >= baseline] -= baseline
    image = Image.fromarray(image)
    
    return image

def BkgSubstraction(image, bkg_cut, show = False):
    """
    We remove the peak from the data and get the baseline. After seeing the 
    baseline shape we perform a 2 order polynomial fit to substract bkg.
    ONLY FOR THE SPECTRUM SHAPE -- not the full image
    """
    #First we show the profile over one axis
    if type(image) == Image.Image:
        data_x = np.sum(image, axis = 0)
    else:
        data_x = image
    x = np.linspace (1, len(data_x), len(data_x)) - np.argmax(data_x)

    #Defining function (2nd order)
    def poly2(x, a, b, c):
        return a*x**2 + b*x + c
    #Auxiliar array where we remove the peak
    data_bg = data_x.astype(dtype= float)
    flat_vals = np.where(data_bg >= min(data_bg) + min(data_bg) * bkg_cut)               #FIXME ---> do it with np.diff or similar
    
    for i in flat_vals: data_bg[i] = np.nan
    #Masiking values != nan so curve_fit can work with them
    mask = np.isfinite(data_bg)
    #Filter x and data_bg using the mask to avoid nan's
    x_filtered = x[mask]
    data_bg_filtered = data_bg[mask]
    #Fit
    opt, cov = curve_fit(poly2, x_filtered, data_bg_filtered)
    #Show
    if show == True:
        plt.figure()
        plt.plot(x, data_x)
        plt.plot(x, poly2(x, *opt))
    #Data wo bkg
    data_substracted = data_x - poly2(x, *opt)
    
    return data_substracted


def ImageCut(geometry, image, vpix, hpix = 0, axis_cut = 'x', offcenter = None, show = False, cmap = 'jet', save = False):
    #geometry = [x, y, width, height]
    if offcenter == None:
        x_offcenter = 0 ; y_offcenter = 0
    else:
        x_offcenter, y_offcenter = offcenter
        
    if axis_cut == 'x':
        x, y, width, height = geometry
        
        left = (width + x_offcenter) / 2 - hpix
        top =  0
        right = (width + x_offcenter) / 2 + hpix
        bottom =  height
         
        cut = image.crop((left, top, right, bottom))
        array_aux = np.array(cut)
        
        imageArray = np.sum(array_aux, axis = 1)
        
        rect = patches.Rectangle((left, bottom),  right - left, top - bottom, 
                                 linewidth=2, edgecolor='r', facecolor='none')
        
    elif axis_cut == 'y':
        x, y, width, height = geometry
        
        left = 0 + hpix + 5
        top =  (height - y_offcenter)/ 2 - vpix
        right = width - hpix +5
        bottom =  (height - y_offcenter )/ 2 + vpix
         
        cut = image.crop((left, top, right, bottom))
        array_aux = np.array(cut)
        
        imageArray = np.sum(array_aux, axis = 0)
        
        rect = patches.Rectangle((left, bottom), right - left, top - bottom, 
                                 linewidth=2, edgecolor='r', facecolor='none')
        
    elif axis_cut == 'square':
        x, y, width, height = geometry
        
        left = (width + 8) / 2 - hpix
        top =  0
        right = (width + 8) / 2 + hpix
        bottom =  height
         
        cut = image.crop((left, top, right, bottom))
        array_aux = np.array(cut)
        
        imageArray = np.sum(array_aux, axis = 1)
        
        rect = patches.Rectangle((left, bottom),  right - left, top - bottom, 
                                 linewidth=2, edgecolor='r', facecolor='none')
    
    if show == True:
        fig, ax = plt.subplots(1)
        a = ax.imshow(image, cmap = cmap)
        ax.add_patch(rect)
        ax.vlines(width/2 + x_offcenter, 182.5, 0)
        ax.tick_params(axis='both', which='major', labelsize = 18)
        fig.colorbar(a)
        if save == True:
            fig.savefig('../tiffs/July24/new set data/highgain-2structures/12x12/100ms/cmap2/ss_single_' + str(i) + '.tif')
    return imageArray, array_aux


def LoadOpticalParams(binning, conf = '_test' ,spath = '../setup_params/datos'):
    path = spath + str(binning) + 'x' + str(binning) + conf + '.csv'
    data = pd.read_csv(path)
    
    return data


def HotPixel12x12(data_image):
    """
    Hot pixel substraction for a 12x12 rebinning.
    Only the image input is needed.
    """
    #FIXME -> create a txt with the pixels positions for each rebinning
    hotpix12x12 = ([9,22], [3,71], [34,2], [23,185], [71,38], [23,185],
                   [154,188], [166,162], [171,109], [167,94], [127,78],
                   [144,20], [71,38], [74,98], [34,117], [56,140])
    
    data_imageNAN = data_image.astype(float)
    #Set hot pixels to nan and calculate the mean without them
    for i in hotpix12x12: data_imageNAN[i[0], i[1]] = np.nan
    mask = np.isfinite(data_imageNAN)
    data_filtered = data_image[mask]
    mean_aux = np.mean(data_filtered)
    #Replace initial values from hotpix to the actual mean
    for i in hotpix12x12: data_image[i[0], i[1]] = mean_aux

    return data_image


def gainFunc(data, qeff=1,geomeff=1,T=1,E=5.5,W=25.7e-6,A=500, exp_time=1, base = 100):
    
    '''
    Calculation of the optical gain given a set of data and the optical 
    parameters and time exposure.
    '''  
    #Remove 100 ADU baseline
    data[data < base] = 0
    data[data >= base] -= base
    #Convert to photons --> 1.32 is the gain from the camera (see manual)
    data = data * 1.32 
    #Correct for the optical parameters
    data = data / qeff /geomeff / T
    total_photons = np.sum(data)
    #Total number of e* produced for the time exposure
    electrons = exp_time * A * E / W
    
    gain = total_photons / electrons
    print('Gain: %.3e \t'%gain)
    
    return gain
