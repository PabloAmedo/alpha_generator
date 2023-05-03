# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:48:14 2023

@author: jacob
"""

from Alpha_track_simulator import*
from optical_gain import*


#==This is the main program==
#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)
#Create a source and a list to store

source=Source(radius=0.1)
track_list=[]

#Creat some alpha tracks and set some parameters
n_tracks=3000;ath_angle=0

exposition_time=n_tracks/source.rate #In seconds


source.produce_alpha(n=n_tracks,store=track_list,ath_in=None,ionization_profile="Bragg")


#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler()
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()


#Create a noise object with a given dark noise
noise=Noise(50)

#Plot the tracks
image2d=Image_2D(track_list=track_list,hist_args={"bins":20})
# #Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object=noise,exposition_time=exposition_time)
image2d.plot_x()





path='C:/Users/jacob/OneDrive - Universidade de Santiago de Compostela/Documentos (Escritorio)/Física/0. Laboratorio DGD/Imágenes Teledyne/'
file='ss_single_1-1_CORTE.tif'
datos='datos.csv'
datos=pd.read_csv(path+datos)

cameraImage=optical_gain(path,file)
cameraImage.plot_data(RCSDA=float(datos['RCSDA']),
                                    Rgem=float(datos['Rgem']),
                                    Rtubo=float(datos['Rtubo']),
                                    pie=100,
                                    cal=float(datos['cal_fab']),
                                    qeff=float(datos['qeff']),
                                    geomeff=float(datos['geomeff']),
                                    T=float(datos['T']))

cameraImage.data_plot=cameraImage.data_plot-np.min(cameraImage.data_plot)

fig,ax1=plt.subplots(ncols=1)
ax1.hist(image2d.electron_cut,bins=100,fill=False,label='Montecarlo')
ax1.set_xlabel('x (cm)')
ax1.set_ylabel('N')

ax1.plot(cameraImage.x,cameraImage.data_plot/np.max(cameraImage.data_plot)*max(image2d.le_hist[0]))

#%%















