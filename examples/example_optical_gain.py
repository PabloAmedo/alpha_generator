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
import matplotlib.patches as patches

print('Running...')
sys.path.append('/..')
start = time.time()
plt.close('all')

save = 0

n_tracks = 100
electrons = int(5.5e6 / 25.7) #214007
ph_target = 2002

red = 100#electrons / ph_target

conf = '{}t_1e{}e'.format(n_tracks, red)
diff = 0.25 
ath_angle = 40 * np.pi / 180
phi_angle = 190 * np.pi / 180
bins = [256,256]
pixel_size = (16e-4 * 4 * 18)

x_range = [-(pixel_size * bins[0])/2, (pixel_size * bins[0])/2]
y_range = [-(pixel_size * bins[1])/2, (pixel_size * bins[1])/2]


#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)            #THIS IS NOT BEING USED AT ALL

#Create a source and a list to store
source=Source(radius = 0.35, red_fact = red)
track_list=[]

exposition_time = n_tracks/source.rate #In seconds
source.produce_alpha(n = n_tracks, store = track_list, ionization_profile = "Bragg", ath_in = ath_angle, phi_in = phi_angle)

#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler(sigma_diff = diff, sigma_PSF = 0)
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()

#Create a noise object with a given dark noise
noise = Noise(110)

#Plot the tracks
image2d = Image_2D(track_list = track_list, hist_args = {"bins": bins, "range": [x_range, y_range]})          #Calcualted to match 12x12 rebinning px size (4.5*12 (um))
# #Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x(variable_cut = 'x', min_var = -0.648, max_var = 0.378)
os.chdir('alpha_generator/')


plt.figure()

plt.imshow(image2d.Hist2D, origin = 'lower', aspect='auto', extent=[x_range[0], x_range[1], y_range[0], y_range[1]], cmap = 'gray')
#circle = patches.Circle((0, 0), 50 * 0.1242, color='red', fill=False, linewidth=2)
#plt.gca().add_patch(circle)

#circle2 = patches.Circle((0, 0), 71 * 0.1242, color='green', fill=False, linewidth=2)
#plt.gca().add_patch(circle2)

plt.tick_params(axis='both', which='major', labelsize = 18)

if save == 1:
    np.savetxt('data_debug/eric/sim' + conf + '.txt', image2d.Hist2D)

end = time.time()
print('Time elapsed:\t', end-start)
#noise es un objeto con dark current = 50 (linea 39)
plt.close('all')

plt.figure()
plt.plot(np.sum(image2d.Hist2D, axis = 0))
plt.xlim([85,135])
