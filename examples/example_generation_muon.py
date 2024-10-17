# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:16:43 2023

@author: diego
"""
import os
os.chdir('../')
from Alpha_track_simulator import *
import scipy.special as spc
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


plt.close('all')
print('Running...\n')

def Momentum(energy, mass):
    p_list = []
    p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)

#INPUTS========================================================================
n_tracks    =   5
P           =   1                                                         #bar
#y0          =   10
#ath0        =   50 * ( np.pi / 180)
E           =   2500                                                          #E = 4000 [MeV] (CR)
dimensions  =   [14,14,14]
mass        =   105.66
gas         =   'ArCF4_99-1'
e_cut       =   10000

px          =   13e-3 
rebin       =   4 
mag         =   11 
px_size     =   (px * rebin * mag) * 1e-1

sigma_diff  =   0.21    #!!!!!!!!!!!!!!!!!! NOT SURE
sigma_PSF   =   0
line        =   False

bins        =   [int(dimensions[0] / px_size), int(dimensions[1] / px_size)]
x_range = [0, (px_size * bins[0])]
y_range = [0, (px_size * bins[1])]

muons       =   [] 

#==============================================================================
#==============================================================================

#Muons generation
muon = muon_generator(energy = E, geometry = dimensions, gas = gas, pressure = P)            #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons, e_cut = e_cut)#, phi_in = 45 * np.pi / 180)#, ath_in =  np.pi/1.5 )     #generate muon's obj stored in muons list

argon = Gas(gas = 'ArCF4', L_drift = abs(dimensions[2] - muons[0].z0), Pressure = P)

#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff = sigma_diff, sigma_PSF=0)
diff_handler.diffuse(muons)
#Filling tracks with e-
for muon in muons: muon.fill(diff = True)

#Noise object
noise=Noise(10)
#Plot the final tracks
image2d=Image_2D(track_list = muons, hist_args={"bins":bins, "range": [x_range, y_range]}, pixel_size = px_size)
#image2d.track_plot()
#image2d.plot_hist(noise_object=noise)
image2d.plot_x()

plt.close('all')

#Getting histogram and edges to plot with white background
hist = image2d.Hist2D_e
hist_x = image2d.x_edges
hist_y = image2d.y_edges

hist[hist == 0] = 'NaN'                                                        #pcolormesh plots nans as white

plt.figure()
plt.pcolormesh(hist_x, hist_y, hist, cmap = 'jet')
plt.colorbar()
plt.xlim([0, dimensions[0]])
plt.ylim([0, dimensions[1]])
plt.xlabel('X (cm)', fontsize = 15)
plt.ylabel('Y (cm)', fontsize = 15)

print('px size:\t({:.3f} x {:.3f}) (cm x cm)'.format(np.mean(np.diff(hist_x)), np.mean(np.diff(hist_y))))

#PLOT =========================================================================
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = []
y = []
z = []
for i in range(len(muons)):
    x = [float(muons[i].x0), float(muons[i].x)]
    y = [float(muons[i].y0), float(muons[i].y)]
    z = [float(muons[i].z0), float(muons[i].z)]

    ax.plot(x, y, z,'-o', color = 'r')
cube_size = dimensions[0]
for i in range(2):
    for j in range(2):
        ax.plot([0, cube_size], [i*cube_size, i*cube_size], [j*cube_size, j*cube_size], color='k')
        ax.plot([i*cube_size, i*cube_size], [0, cube_size], [j*cube_size, j*cube_size], color='k')
        ax.plot([i*cube_size, i*cube_size], [j*cube_size, j*cube_size], [0, cube_size], color='k')
        
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')

ax.set_xlim([0 - 3 , dimensions[0] + 3])
ax.set_ylim([0 - 3 , dimensions[1] + 3])
ax.set_zlim([0 - 3, dimensions[2] + 3])
ax.grid(False)
plt.show()