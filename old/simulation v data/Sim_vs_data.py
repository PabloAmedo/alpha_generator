# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 12:08:46 2023

@author: diego
"""
import os
print('Running...')
os.chdir('../')
print(os.getcwd())
import numpy as np
from Alpha_track_simulator import*
from PIL import Image
import tifffile
from Alpha_track_simulator import*
from optical_gain import*


"""
Testing Simulation and comparing results with measured data
"""


n_tracks=1              #tracks to generate --> exposure time is obtained from here
reduc_fact=1          #factor to reduce the number of generated e-
rebin='12x12'           #rebin performed in the real image used
scale=1.0               #scale to match real vs sim image (only with higher rebin ?????)



#%% ============ GENERATION ============

#Generate source

source=Source(radius=0.1)

#Generate alpha track

tracks=[]
source.produce_alpha(n=n_tracks, store=tracks, ionization_profile='Bragg',n_electrons=int(214007/reduc_fact))             #TRACK PRODUCTION --- change n to number of tracks

exposition_time=n_tracks/source.rate


#%% ============ DIFFUSION ============

diff=Diffusion_handler(sigma_diff=0.25,sigma_PSF=0)                             #Data from Jacobo¡s TFG
diff.diffuse(tracks)

# We generate the primary electrons + diffused ones

for i in tracks: i.fill()


#%% ============ 1st PLOT ============

#Plot the tracks
image2d = Image_2D(track_list=tracks,hist_args={"bins":500})   #default was 100

noise=Noise(50)                                                                 #Default data in the original example 
# #Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object=noise,exposition_time=exposition_time)
image2d.plot_x()


#%% ============ NOISE + SAVE + LOAD ============
#Adding some noise and SAVING image

IMAGENtiff=noise.add_noise(10, image2d.Hist2D)  #default was 10
path='alpha_generator/tiffs/'
image=Image.fromarray(IMAGENtiff)
image.save(path+'sim_image.tiff')

#Loading real image

file='ss_single_3_1s.tiff' #'ss_single_1-1_CORTE.tif'
datos='../setup_params/datos' + rebin +'.csv'                                                          # Use datos_custom to test new set ups
datos=pd.read_csv(path+datos)

cameraImage=data_image(path + file)
cameraImage.plot_data(RCSDA=float(datos['RCSDA']),
                      Rgem=float(datos['Rgem']),
                      Rtubo=float(datos['Rtubo']),
                      cal=float(datos['cal_fab'])
                      )


#%% ============ DATA MANIPULATION ============

# Obtain the gain
cameraImage.gain(qeff=float(datos['qeff']),geomeff=float(datos['geomeff']),T=float(datos['T']), reduc_fact=reduc_fact, exp_time=exposition_time)


# PHOTONS IN X-AXIS

# Filter the x data taken from the camera to match the alpha range

r=4.8

camera_x=[x for x in cameraImage.x if x >= -r and x <= r]                      # we select the values from the REAL image where the alpha could be spotted
                                                                               # cameraImage.x are the values in cm for the pixels/bins ?? of the camera

positions=[(i,x) for (i,x) in enumerate(cameraImage.x) if x >= -r and x <= r]  # enumerate(camera_x)]    #is the same, right??? 
positions=[x[0] for x in positions]

# Take the accumulated numb of photons in x-axis (stored in .data_to_plot_x)

photons_x=[cameraImage.data_to_plot_x[i] for i,_ in enumerate(cameraImage.data_to_plot_x) if i in positions]
#photons_x=photons_x/max(image2d.le_hist[0])

#%% ============ REBIN(?) & BG SUBSTRACTION ============

# REBIN (???)

# First we are going to set a rebinning process
"""
step=1
# Uncomment below section when a rebinning is performed --->  need to change definition of photons_x

# Check to improve this in a future
aux_photons=[]
aux_x=[]

for i in range(0,len(photons_x)-step,step):
    
    aux_photons.append(np.sum(photons_x[i:i+step]))
    aux_x.append(np.sum(camera_x[i:i+step])/step)

photons_x=np.copy(aux_photons)
camera_x=np.copy(aux_x)

"""
# BG 'SUBSTRACTION'

# Removing bg linear contribution. First we model the linear function
# I think this part could be improved!!

b= ( (np.mean(photons_x[0])-np.mean(photons_x[-1]))/(np.mean(camera_x[0])-np.mean(camera_x[-1])) )

a=( np.mean(photons_x[0])-b*np.mean(camera_x[0]) )

for i in range(len(photons_x)):
    
    photons_x[i]=photons_x[i]-(a+b*camera_x[i])
    
    if photons_x[i]<0:
        photons_x[i]=0
        
#%% ============ PLOTS ============


fontsize=20
fig,ax1=plt.subplots(ncols=1)


ax1.set_title('Data vs Montecarlo',fontsize=fontsize)
ax1.hist(image2d.electron_cut,bins=100,fill=False,label='Montecarlo')
ax1.set_xlabel('x (cm)',fontsize=fontsize)
ax1.set_ylabel('N',fontsize=fontsize)
ax1.grid(True)



ax1.plot(camera_x,(photons_x/np.max(photons_x))*max(image2d.le_hist[0]),color='blue',label='Data')
ax1.legend()

plt.figure()
plt.imshow(cameraImage.data, cmap='gray')
plt.title('REAL IMAGE -- {} s'.format(exposition_time))


"""     Some interesting plots that can be used to test things


plt.figure()
plt.imshow(IMAGENtiff, cmap='gray')


plt.figure()

plt.plot(camera_x,(photons_x/np.max(photons_x))*max(image2d.le_hist[0]),color='blue',label='Data')

"""

# FINAL IMAGE GENERATION

x_size=cameraImage.data.shape[0]
y_size=cameraImage.data.shape[1]
media=(np.mean(IMAGENtiff[0])+np.mean(IMAGENtiff[1]))/2
#finalImage=np.random.uniform(low=0.98, high=1.01 ,size=(int(1.5*x_size),int(1.5*y_size)))*media
finalImage=np.zeros((int(scale*x_size),int(scale*y_size)))
fImage=noise.add_noise(10, finalImage)


fImage[int(scale*x_size/2-50):int(scale*x_size/2+50),int(scale*y_size/2-50):int(scale*y_size/2+50)]=IMAGENtiff
fscale=np.around(1/scale,2)

plt.figure()


plt.imshow(fImage, cmap='gray')
plt.title('SIMULATED IMAGE -- {} tracks // scale=x{}'.format(n_tracks,fscale))


