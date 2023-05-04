# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:48:14 2023

@author: jacob
"""

from Alpha_track_simulator import*
from optical_gain import*

plt.close('all')

#==This is the main program==
#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)
#Create a source and a list to store

source=Source(radius=0.1)
track_list=[]

#Creat some alpha tracks and set some parameters
n_tracks=20000;ath_angle=0

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

# cameraImage.data_plot=cameraImage.data_plot-np.min(cameraImage.data_plot)
r=4.8
camera_x = [x for x in cameraImage.x if x <= r and x>= -r]                      # Posiciones entre -5 y 5 cm
positions = [(i,x) for i,x in enumerate(cameraImage.x) if x <= r and x>-r]      # Paso intermedio
positions = [x[0] for x in positions]                                           # Índices de las posiciones

camera_data = [cameraImage.data_plot[i] for i, _ in enumerate(cameraImage.data_plot) if i in positions]

paso=30
new_vector=[]
new_vector_x=[]
for i in range(0,len(camera_data)-paso,paso):
    new_vector.append(np.sum(camera_data[i:i+paso]))
    new_vector_x.append(np.sum(camera_x[i:i+paso])/paso)

camera_data=np.copy(new_vector)
camera_x=np.copy(new_vector_x)

b = (np.mean(camera_data[0])-camera_data[-1])/(camera_x[0]-camera_x[-1])*1.6
a = camera_data[0] - b*camera_x[0]

for i in range(len(camera_data)):
    camera_data[i] = camera_data[i] - (a + b*camera_x[i])
    if camera_data[i]<0:
        camera_data[i]=0


fontsize=20
fig,ax1=plt.subplots(ncols=1)
ax1.set_title('Data vs Montecarlo',fontsize=fontsize)
ax1.hist(image2d.electron_cut,bins=100,fill=False,label='Montecarlo')
ax1.set_xlabel('x (cm)',fontsize=fontsize)
ax1.set_ylabel('N',fontsize=fontsize)
ax1.grid(True)

ax1.plot(camera_x,camera_data/np.max(camera_data)*max(image2d.le_hist[0]),color='blue',label='Data')
ax1.legend()







