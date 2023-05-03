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

class optical_gain:
        
    def __init__(self, path, file):
         
        self.image = Image.open(path+file)
        self.data  = np.array(self.image)
        self.data_plot = np.zeros(np.shape(self.data)[1])
        
    def plot_data(self,RCSDA=1,Rgem=1,Rtubo=1,pie=100,cal=1.,qeff=1,geomeff=1,T=1):
        
        # Quitamos el ruido electrónico
        for i in range(np.shape(self.data)[0]):
            for j in range(np.shape(self.data)[1]):
                if self.data[i,j] < pie:
                    self.data[i,j] = 0
                else:
                    self.data[i,j] = self.data[i,j]-pie
                    
        if qeff==1:
            self.data = self.data
        else:
            self.data = self.data/qeff
    
        if T==1:
            self.data = self.data
        else:
            self.data = self.data/T

        if geomeff==1:
            self.data = self.data
        else:
            self.data = self.data/geomeff
    
        # Multiplicamos por el factor de ganancia de la cámara
        self.data = self.data*1.32

        for i in range(np.shape(self.data)[0]):
            self.data_plot = self.data_plot + self.data[i,:]

        # Centramos el eje x
        self.x = np.arange(0,len(self.data_plot),1) - (np.where(self.data_plot==max(self.data_plot)))[0][0]
        self.x = self.x/cal

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
        plt.plot(self.x, self.data_plot, '-', color='black')
        
        plt.legend(loc='best')
        
      
    def gain(self,data_plot):
        
        self.total_photons = np.sum(data_plot)
        print('El número total de fotones es de %f' %np.sum(self.total_photons))

        E = 5.5e6 ; W = 26.7 ; A = 500
        self.electrons = 30*A*E/W
        self.gain = self.total_photons/self.electrons
        print('Eficiencia fotones/e_primario: %.3f'%self.gain)
        return(self.gain)
    
