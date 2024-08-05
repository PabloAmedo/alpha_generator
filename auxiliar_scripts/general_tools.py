# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:32:00 2024

@author: usuario
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


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