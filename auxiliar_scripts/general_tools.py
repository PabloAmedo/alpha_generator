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

def cmap_change(image, cmap, pct_max = 0.5):
    
    img_array = np.array(image)
    
    fig = plt.figure()
    plt.title('{}'.format(cmap), fontsize = 25)
    img = plt.imshow(image, cmap = cmap, vmax = np.max(img_array) * pct_max)
    plt.colorbar(img)
    plt.tick_params(axis='both', which='major', labelsize = 18)
    plt.xlabel('X (px)', fontsize = 20)
    plt.ylabel('Y (px)', fontsize = 20)
    fig.tight_layout()
    

def posicion_segundo_maximo(arr):
    # Verificar si hay al menos dos elementos únicos
    if len(arr) < 2:
        return None  # No hay segundo máximo

    # Inicializar el máximo y segundo máximo
    maximo = segundo_maximo = float('-inf')
    pos_maximo = pos_segundo_maximo = -1

    # Encontrar el máximo
    for i, num in enumerate(arr):
        if num > maximo:
            maximo = num
            pos_maximo = i

    # Encontrar el segundo máximo
    for i, num in enumerate(arr):
        if num > segundo_maximo and num < maximo:
            segundo_maximo = num
            pos_segundo_maximo = i

    if segundo_maximo == float('-inf'):
        return None  # No hay segundo máximo (por ejemplo, si todos los elementos son iguales)

    return pos_segundo_maximo