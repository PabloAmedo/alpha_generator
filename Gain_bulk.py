# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 13:46:30 2023

@author: diego
"""

import os
from Alpha_track_simulator import*
import tifffile
from optical_gain import*


rebin='12x12'                                                                  #according to image readout
exposition_time=10 #s                                                          #according to image readout
path='C:/Users/usuario/Desktop/PRACTICAS - IGFAE/IGFAE - GASEOUS DETECTORS/AnalysisFramework/ImageAnalysis/alpha_generator/' #main directory
folder='tiffs/Oct_run/10 seconds/1_st_GM_on/'                                                #images location


datos='datos'+rebin+'.csv'                                                     #set-up data readout
datos=pd.read_csv(path+datos)

tiff_files=[f for f in os.listdir(path+folder) if f.endswith('.tiff')]         #select speciffic format for images


for file in tiff_files:
    
    cameraImage=data_image(path+folder,file)
    print('file: ', file)
    cameraImage.gain(qeff=float(datos['qeff']),geomeff=float(datos['geomeff']),T=float(datos['T']), reduc_fact=1, exp_time=exposition_time)
    print('\n')