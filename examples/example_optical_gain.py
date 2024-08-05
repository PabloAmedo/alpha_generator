# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:48:14 2023

@author: jacob
"""
import os
os.chdir('../')

from PIL import Image
import tifffile
from Alpha_track_simulator import*
from optical_gain import*
import sys
import time

print('Running...')
sys.path.append('/..')
start = time.time()
plt.close('all')


n_tracks = 350
red = 1000
conf = '{}t_1e{}e'.format(n_tracks, red)
diff = 0.27
ath_angle = 0

#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)            #THIS IS NOT BEING USED AT ALL

#Create a source and a list to store
source=Source(radius = 0.35, red_fact = red)
track_list=[]

exposition_time = n_tracks/source.rate #In seconds
source.produce_alpha(n = n_tracks, store = track_list, ath_in = None, ionization_profile = "Bragg")

#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler(sigma_diff = diff,sigma_PSF=0)
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()

#Create a noise object with a given dark noise
noise = Noise(50)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins":100})          #Calcualted to match 12x12 rebinning px size (4.5*12 (um))
# #Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x(variable_cut = 'x', min_var = -0.648, max_var = 0.378)
os.chdir('alpha_generator/')
plt.figure()
plt.imshow(image2d.Hist2D)

np.savetxt('data_debug/diff_2-7mm/sim' + conf + '.txt', image2d.Hist2D)

end = time.time()
print('Time elapsed:\t', end-start)
#noise es un objeto con dark current = 50 (linea 39)
"""
IMAGENtiff = noise.add_noise(10,image2d.Hist2D)
image = Image.fromarray(IMAGENtiff)
image.save('simulated_image_to_analyze.tiff')
# tifffile.imwrite(path+'simulated_image_to_analyze.tif', IMAGENtiff)
#==============================================================================
#REAL DATA LOADING

print(os.getcwd())
file ='tiffs/July/c-1136_r-128_acground_aa1950_fc2450_fa3480/50ms/ss_single_1_cut3.tif' #'ss_single_1-1_CORTE.tif'
datos ='setup_params/datos12x12_test.csv'
datos = pd.read_csv(datos)

cameraImage = data_image(file)
cameraImage.plot_data(RCSDA = float(datos['RCSDA']),
                      Rgem = float(datos['Rgem']),
                      Rtubo = float(datos['Rtubo']),
                      cal = float(datos['cal_fab']),
                      )

# La ganancia da sobre 70, pero no estamos en la imagen entera
cameraImage.gain(qeff = float(datos['qeff']),geomeff = float(datos['geomeff']),T = float(datos['T']))

#cameraImage.data_plot=cameraImage.data_plot-np.min(cameraImage.data_plot)
r = 5
camera_x = [x for x in cameraImage.x if x <= r and x>= -r]                      # Posiciones entre -r y r cm
positions = [(i,x) for i,x in enumerate(cameraImage.x) if x <= r and x>-r]      # Paso intermedio
positions = [x[0] for x in positions]       

#Índices de las posiciones
camera_data = [cameraImage.data_to_plot_x[i] for i, _ in enumerate(cameraImage.data_to_plot_x) if i in positions]
#camera_data = camera_data / max(image2d.le_hist[0])                           #only uncomment if we are performing rebinning
"""
# =============================================================================
# =============================================================================
# =============================================================================
"""
#Is this some kind of re-binning ??????
    
paso = 1                                        
new_vector = []
new_vector_x = []

for i in range(0, len(camera_data)-paso, paso):
    new_vector.append(np.sum(camera_data[i:i + paso]))
    new_vector_x.append(np.sum(camera_x[i:i+paso]) / paso)

camera_data = np.copy(new_vector)
camera_x = np.copy(new_vector_x)
"""
# =============================================================================
# =============================================================================
# =============================================================================
"""
# SUBSTRACTING LINEAR CONTRIBUTION (BG ???)

b = ( (np.mean(camera_data[0]) - np.mean(camera_data[-1])) / (np.mean(camera_x[0]) - np.mean(camera_x[-1])) ) #*2.5   #calculo dunha pendente
a = ( np.mean(camera_data[0]) - b * np.mean(camera_x[0]) )  #* 1.2                # ordeada no orixe

#b = ( camera_data[0] - camera_data[-1] ) / (camera_x[0] - camera_x[-1] ) * 2
#a =( camera_data[0] - b*camera_x[0] ) + 150

#camera_data=camera_data[:-1]

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
ax1.plot(camera_x,(camera_data/np.max(camera_data))*max(image2d.le_hist[0]),color='blue',label='Data')
ax1.legend()

end = time.time()
print('Time elapsed:\t', end-start)
"""