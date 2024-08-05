# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:22:42 2024

@author: diego
"""
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import pandas as pd
import os

#FUNCTIONS ---> This can be merged into a class OpticalGain
class DataImage:
    
    def __init__(self, file):#DEFINE HERE --> Load image, maybe setup settings too??
        self.image = Image.open(file)
        self.data = np.array(self.image)
        
        
    def ImageCut(self, geometry, vpix, hpix, offcenter = [], axis_cut = 'y', show = False, path_to_save = None, cmap = None, vmax = None ):
        #geometry = [x, y, width, height]
        x_offcenter, y_offcenter = offcenter
        x, y, width, height = geometry
        if axis_cut == 'x':
        
            left = (width + x_offcenter) / 2 - hpix
            top =  0
            right = (width + x_offcenter) / 2 + hpix
            bottom =  height
             
            cut = self.image.crop((left, top, right, bottom))
            self.data_cut = np.array(cut)
            
            self.data_sum = np.sum(self.data_cut, axis = 1)
            
            rect = patches.Rectangle((left, bottom),  right - left, top - bottom, 
                                     linewidth=2, edgecolor='r', facecolor='none')
            
        elif axis_cut == 'y':
            
            left = 0 + hpix + x_offcenter
            top =  (height - y_offcenter)/ 2 - vpix
            right = width - hpix + x_offcenter
            bottom =  (height - y_offcenter )/ 2 + vpix
             
            cut = image.crop((left, top, right, bottom))
            self.data_cut = np.array(cut)
            
            self.data_sum = np.sum(self.data_cut, axis = 0)
            
            rect = patches.Rectangle((left, bottom), right - left, top - bottom, 
                                     linewidth=2, edgecolor='r', facecolor='none')
            
        elif axis_cut == 'square':
            
            left = (width + x_offcenter) / 2 - hpix
            top =  0 - y_offcenter
            right = (width + x_offcenter) / 2 + hpix
            bottom =  height + y_offcenter
             
            cut = image.crop((left, top, right, bottom))
            self.data_cut = np.array(cut)
            
            self.data_sum = np.sum(self.data_cut, axis = 0)
            
            rect = patches.Rectangle((left, bottom),  right - left, top - bottom, 
                                     linewidth=2, edgecolor='r', facecolor='none')
        
        if show == True:
            fig, ax = plt.subplots(1)
            img = ax.imshow(self.image, cmap = cmap, vmax = vmax)
            #ax.add_patch(rect)
            ax.tick_params(axis='both', which='major', labelsize = 18)
            fig.colorbar(img)
            if save == True:
                fig.savefig(path_to_save + '.tif')
                
        return self.sum, self.cut_data
    
    def LoadSimulation(self, file, path = None): #LOAD SIMULATION!!!!!!!
        if path == None:
            self.simulated_data = np.loadtxt(file)
        else:
            self.simulated_data = np.loadtxt(path + file)
        self.simulated_cut = self.simulated_data[:][26:76]
        self.simulated_sum = np.sum(self.simulated_cut, axis = 0)
    
    def LoadOpticalParams(self, binning, conf = '_test' ,spath = '../setup_params/datos'):
        """
        This method loads some optical parameters related to the setup 
        configuration from a .txt file into a pandas dataframe.
        """
        path = spath + str(binning) + 'x' + str(binning) + conf + '.csv'
        self.optdata = pd.read_csv(path)
        
        return self.optdata
    
    def BaselineSubstraction(self, baseline = 100):                     #This is good, maybe no need a func
        """
        Substracts the baseline value for everypixel. The default value is an 
        estimation, according to previous analysis the average value 
        is 108.9 e/px
        """
        self.data[self.data < baseline] = 0
        self.data[self.data >= baseline] -= baseline
        image = Image.fromarray(self.data)
        
        return image
    
    def Gain(self):
        nedata = sum(self.data_sum) / float(self.optdata['qeff']) / float(self.optdata['geomeff']) / float(self.optdata['T'])
        nesim = sum(h2sum)
        gain = nedata / nesim