# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 13:46:30 2023

@author: diego
"""

import os
from Alpha_track_simulator import*
import tifffile
from optical_gain import*


rebin='2x2'
exposition_time=1
path='C:/Users/diego/Desktop/PRACTICAS-IGFAE/AnalysisFramework/Simulation/alpha_generator/'
folder='tiffs/'


datos='datos'+rebin+'.csv'
datos=pd.read_csv(path+datos)

tiff_files=[f for f in os.listdir(path+folder) if f.endswith('.tiff')]


for file in tiff_files:
    
    cameraImage=data_image(path+folder,file)

    cameraImage.gain(qeff=float(datos['qeff']),geomeff=float(datos['geomeff']),T=float(datos['T']), reduc_fact=1, exp_time=exposition_time)