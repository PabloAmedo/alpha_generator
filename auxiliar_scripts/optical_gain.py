# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 14:15:15 2023

@author: jacob
"""

"""---------------------------------------------------------------"""
"""---------------------------------------------------------------"""
"""Load an image and get its optcial gain, plot, and total photons"""
"""---------------------------------------------------------------"""
"""---------------------------------------------------------------"""

from PIL import Image
from Alpha_track_simulator import Source, Alphas_tracks
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class data_image():
    
    '''This class loads an image and plot it'''
    
    def __init__(self,file, path=None,pie = 100):
        if path != None:
            self.image = Image.open(path+file)
        else:
            self.image = Image.open(file)
        self.data  = np.array(self.image)
        self.data_plot = np.zeros(np.shape(self.data)[1])
        
        # Quitamos el ruido electrónico
        # for i in range(np.shape(self.data)[0]):
        #     for j in range(np.shape(self.data)[1]):
        #         if self.data[i,j] < pie:
        #             self.data[i,j] = 0
        #         else:
        #             self.data[i,j] = self.data[i,j] - pie
        #-------------------------------------------------
        self.data[self.data < pie] = 0
        self.data[self.data >= pie] -= pie     
        #-------------------------------------------------                    
        '''Quitamos el ruido electrónico en la propia inicialización
        de la imagen, en caso de no querer hacerlo simplemente pie=0
        '''

    def plot_data(self,RCSDA=1.,Rgem=1.,Rtubo=1.,cal=1.):                       #Params inicializados, seran leidos desde archivo
        
        '''
        Esta función es ÚNICAMENTE para plotear los datos en x
        '''
                
        # for i in range(np.shape(self.data)[0]):
        #     self.data_plot = self.data_plot + self.data[i,:]
        #-------------------------------------------------
        self.data_to_plot_x = np.sum(self.data, axis = 0)                          #Valores acumulados para cada pix en x
        #-------------------------------------------------

        # Centramos el eje x
        self.x = np.arange(0,len(self.data_to_plot_x),1) - (np.where(self.data_to_plot_x == max(self.data_to_plot_x)))[0][0]
        self.x = self.x / cal


        plt.figure(figsize=(10,6),dpi=120)
        plt.title('Número de fotones totales en función de x')
        plt.grid(True)

        # Centro x=0
        plt.axvline(x=0, color='black', linestyle='--')
    
        # Límites de el rango de la partícula alpha
        plt.axvline(x=RCSDA, color='green',linestyle='-',label=u'$R_{CSDA}$')
        plt.axvline(x=-RCSDA,color='green',linestyle='-')
    
        # Límites del fatGEM
        plt.axvline(x=Rgem, color='red',linestyle='-',label=u'$fatGEM$')
        plt.axvline(x=-Rgem,color='red',linestyle='-')
    
        # Límites del tubo    
        plt.axvline(x=Rtubo, color='blue',linestyle='-',label=u'$Tubo$')
        plt.axvline(x=-Rtubo,color='blue',linestyle='-')

        if cal>1:
           plt.xticks(np.arange(int(min(self.x)),int(max(self.x)),1))
           
        plt.ylabel('photons',size=10)
        plt.xlabel('x(cm)',size=10)
        plt.plot(self.x, self.data_to_plot_x, '-', color='black')
        
        plt.legend(loc='best')
    
      
    def gain(self,qeff=1,geomeff=1,T=1,E=5.5,W=25.7e-6,A=500, exp_time=1):
    
        '''
        Calculation of the optical gain given a set of data and the optical 
        parameters and time exposure.
        '''           
        self.data = self.data / qeff / geomeff / T
        #ADU to ph
        self.data = self.data * 1.32                                           # 1.32 is the camera gain
        self.total_photons = np.sum(self.data)
        #Number of e- produced for the given exp time
        self.electrons = exp_time * A * E / W
        
        self.gain = self.total_photons/self.electrons
        print('Gain: %.3e \t'%self.gain)
        
        return self.gain

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
        
