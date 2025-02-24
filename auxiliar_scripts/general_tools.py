# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:32:00 2024

@author: usuario
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import PIL


#DEFINING SOME CONVERSION CONSTANTS
deg_to_rad =  np.pi / 180


## FLUCTUATIONS IN e-'s READOUT
def statistics_avalanche(e_per_cluster):
    nClusters = len(e_per_cluster)
    detectedElectrons   = []
    
    for i in range(nClusters):
        acum = 0
        for _ in range(int(e_per_cluster[i])):
            acum = acum + np.random.exponential(scale = 1 )
        detectedElectrons.append(acum) #FIXME

    return np.array(detectedElectrons)



#GENERAL FUNCTIONS
def generate_gif(folder, time_step = 0.25, path = None, ext = '.tiff', nametosave = None, save = False):
    
    if path == None:
        im_path = ''
    else:
        im_path = path

    tiff_files = [f for f in os.listdir(im_path + folder) if f.endswith(ext)]
    
    def number_sort(file_name):
        match = re.search(r'(\d+)', file_name)
        return int(match.group(1)) if match else 0

    sorted_tiff_files = sorted(tiff_files, key = number_sort)

    images = []
    for tiff in sorted_tiff_files:
        images.append(imageio.imread(im_path + folder + tiff))
    if save == True:
        imageio.mimsave(nametosave + '.gif', images, duration = time_step)
    
    return images

def cmap_change(image, cmap, pct_max = 1, title = None):
    img_array = np.array(image)
    
    fig = plt.figure()
    if title != None:
        plt.title('{}'.format(title), fontsize = 18)
    img = plt.imshow(image, cmap = cmap, vmax = np.max(img_array) * pct_max)
    plt.colorbar(img)
    plt.tick_params(axis='both', which='major', labelsize = 13)
    plt.xlabel('X (px)', fontsize = 15)
    plt.ylabel('Y (px)', fontsize = 15)
    fig.tight_layout()