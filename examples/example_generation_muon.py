# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:16:43 2023

@author: diego
"""
import os
os.chdir('../')
from Alpha_track_simulator import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import time
import matplotlib

plt.close('all')
print('Running...\n')
t0 = time.time()
def Momentum(energy, mass):
    p_list = []
    p_list.append(np.sqrt((energy + mass)**2 - mass**2))
        
    return np.array(p_list)

#INPUTS========================================================================
n_tracks    =   3
P           =   1.5                                                            #bar
#y0          =   10
#ath0        =   50 * ( np.pi / 180)
E           =   3500                                                           #E = 4000 [MeV] (CR)
dimensions  =   [16,16,19]
mass        =   105.66
gas         =   'ArCF4_99-1'
e_cut       =   10000
gain        =   2000
#Camera settings
px          =   13e-3 
rebin       =   4 
mag         =   11
sigma_bkg   =   0.3 
counts_bkg  =   1.5

px_size     =   (px * rebin * mag) * 1e-1

sigma_diff  =   0.21    #can be obtained from diff_handler itself
sigma_PSF   =   0
line        =   False

bins        =   [int(dimensions[0] / px_size), int(dimensions[1] / px_size)]
x_range = [0, (px_size * bins[0])]
y_range = [0, (px_size * bins[1])]

muons       =   []
p3D         =   1

#==============================================================================
#==============================================================================

#Muons generation
muon = muon_generator(energy = E, geometry = dimensions, gas = gas, pressure = P)            #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons, e_cut = e_cut)#,
                  #position_in = [25, 10, 12.5])#, ath_in = 90 * np.pi / 180)#, ath_in = 100 * np.pi / 180 )     #generate muon's obj stored in muons list

argon = Gas(gas = 'ArCF4', L_drift = abs(dimensions[2] - muons[0].z0), Pressure = P)

#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff = sigma_diff, sigma_PSF=0)
diff_handler.diffuse(muons)
#Filling tracks with e-
for muon in muons: muon.fill(diff = True)

#Noise object
noise=Noise(counts_bkg)
#Plot the final tracks
image2d = Image_2D(track_list = muons, hist_args={"bins":bins, "range": [x_range, y_range]},
                   gain = gain, QE = 0.95, GE = 5.09e-4)
#image2d.track_plot()
#image2d.plot_hist(noise_object=noise)
#image2d.plot_x()
img_fromCamera = noise.add_noise(0.1, image2d.Hist2D, sigma = sigma_bkg) 

#Getting histogram and edges to plot with white background
hist = image2d.Hist2D
hist_x = image2d.x_edges
hist_y = image2d.y_edges

hist[hist == 0] = 'NaN'    
                                                   #pcolormesh plots nans as white
plt.figure()
plt.imshow(img_fromCamera, cmap = 'gray')
plt.colorbar()
"""
plt.figure()
plt.pcolormesh(hist_x, hist_y, hist, cmap = 'gray')
plt.colorbar()
plt.xlim([0, dimensions[0]])
plt.ylim([0, dimensions[1]])
plt.xlabel('X (cm)', fontsize = 15)
plt.ylabel('Y (cm)', fontsize = 15)
"""
print('px size:\t({:.3f} x {:.3f}) (cm x cm)'.format(np.mean(np.diff(hist_x)), np.mean(np.diff(hist_y))))

#PLOT =========================================================================
if p3D == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = []
    y = []
    z = []
    colors = colors = list(matplotlib.colors.XKCD_COLORS)
    for i in range(len(muons)):
        x = muons[i].electron_positions_diff[:,0]
        y = muons[i].electron_positions_diff[:,1]
        z = muons[i].electron_positions_diff[:,2]
    
        ax.plot(x, y, z,'.', color = colors[i])
    cube_size = dimensions[0]
    for i in range(2):
        for j in range(2):
            ax.plot([0, dimensions[0]], [i*dimensions[1], i*dimensions[1]], [j*dimensions[2], j*dimensions[2]], color='k')
            ax.plot([i*dimensions[0], i*dimensions[0]], [0, dimensions[1]], [j*dimensions[2], j*dimensions[2]], color='k')
            ax.plot([i*dimensions[0], i*dimensions[0]], [j*dimensions[1], j*dimensions[1]], [0, dimensions[2]], color='k')
            
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')
    max_dim = max(dimensions)
    ax.set_xlim([0 - 3 , max_dim + 3])
    ax.set_ylim([0 - 3 , max_dim + 3])
    ax.set_zlim([0 - 3, max_dim + 3])
    ax.grid(False)
    plt.show()

tf = time.time()

print('Time elapsed:\t{} s'.format(tf - t0))