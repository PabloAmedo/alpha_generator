# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 12:07:48 2024

@author: usuario
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
import warnings
from general_tools import *

warnings.filterwarnings("ignore", category=DeprecationWarning)

print('Running...')
sys.path.append('/..')
start = time.time()
plt.close('all')

#INPUTS =======================================================================

n_tracks = 1
red = 10

conf = '{}t_1e{}e'.format(n_tracks, red)
diff = 0.23 /2
ath_angle = 37 * np.pi / 180
phi_angle = 42 * np.pi / 180
bins = [224,183]
pixel_size = (4.5e-4 * 12 * 23)

x_range = [-(pixel_size * bins[0])/2, (pixel_size * bins[0])/2]
y_range = [-(pixel_size * bins[1])/2, (pixel_size * bins[1])/2]

x_offcenter_px = -4
y_offcenter_px = -9

x_offcenter = x_offcenter_px * pixel_size
y_offcenter = y_offcenter_px * pixel_size

#Load data --------------------------------------------------------------------
data = np.loadtxt('Gain analysis/ImageJ_results/Result of ss_single_13 12px.csv', delimiter = ',', skiprows = 1)

data_x, data_y, data_v = np.split(data, 3, axis=1)

data_x = np.array([float(x) for x in data_x])
data_y = np.array([bins[1] - float(y) for y in data_y])
data_v = np.array([float(v) for v in data_v])

#SIMUL GENERATION =============================================================
source=Source(radius = 0.35, red_fact = red)
track_list=[]

exposition_time = n_tracks/source.rate #In seconds
source.produce_alpha(n = n_tracks, store = track_list, phi_in = phi_angle, ath_in = ath_angle, ionization_profile = "Bragg")

diff_handler=Diffusion_handler(sigma_diff = diff, sigma_PSF = 0)
diff_handler.diffuse(track_list)

for i in track_list: i.fill()

noise = Noise(110)

#Plot the tracks --------------------------------------------------------------
image2d = Image_2D(track_list = track_list, hist_args = {"bins": bins, "range": [x_range, y_range]})          #Calcualted to match 12x12 rebinning px size (4.5*12 (um))

#image2d.track_plot()
#image2d.plot_hist(noise_object = noise, exposition_time = exposition_time)
image2d.plot_x(variable_cut = 'x', min_var = -0.648, max_var = 0.378)
os.chdir('alpha_generator/')

img = image2d.Hist2D

plt.figure()
plot = plt.imshow(img,  cmap = 'gray')
plt.plot(data_x + x_offcenter_px , data_y + y_offcenter_px, 's', alpha = 0.5, zorder = 3)

plt.figure()
plot = plt.imshow(img, origin = 'lower', aspect='auto', extent=[x_range[0], x_range[1], y_range[0], y_range[1]], cmap = 'gray')
circle = patches.Circle((0, 0), 50 * pixel_size, color='red', fill=False, linewidth=2)
plt.gca().add_patch(circle)
circle2 = patches.Circle((0, 0), 71 * pixel_size, color = 'green', fill = False, linewidth = 2)
plt.gca().add_patch(circle2)
plt.tick_params(axis='both', which='major', labelsize = 18)
plt.colorbar(plot)

sim_corrected_x = data_x * pixel_size - bins[0]* pixel_size/2 + x_offcenter
sim_corrected_y = data_y * pixel_size - bins[1]* pixel_size/2 + y_offcenter

plt.figure()
plot = plt.imshow(img, origin = 'lower', aspect='auto', extent = [x_range[0], x_range[1], y_range[0], y_range[1]], cmap = 'gray')
plt.plot(sim_corrected_x , sim_corrected_y, 's', alpha = 0.5, zorder = 3)
plt.colorbar(plot)

#ANALYSIS ====================================================================

cut_coord = np.vstack((data_x + x_offcenter_px, data_y + y_offcenter_px)).T
cut_coord = cut_coord.astype(int)

x_min = int(min(data_x)) ; x_max = int(max(data_x))
y_min = int(min(data_y)) ; y_max = int(max(data_y))


#size_x = x_max - x_min ; size_y = y_max - y_min
size_x = 224 ; size_y = 183


data_img = np.zeros((size_x, size_y))

for x, y, adu in zip(data_x, data_y, data_v):
    if 0 <= x < size_x and 0 <= y < size_y:  
        data_img[int(x), int(y)] = adu

sim_img = np.zeros((size_x, size_y))

for x, y in zip(data_x + x_offcenter_px, data_y + y_offcenter_px):
    if 0 <= x < size_x and 0 <= y < size_y:  
        sim_img[int(x), int(y)] = img[int(y), int(x)]

#sim_img[:, 1] = sim_img[:, 1][::-1]

data_img = data_img[x_min:x_max, y_min:y_max]
sim_img = sim_img[x_min:x_max, y_min:y_max]

plt.figure()
plt.imshow(data_img)

plt.figure()
plt.imshow(sim_img)

data_profile = np.sum(data_img, axis = 0)
data_sum = sum(data_profile)

sim_profile = np.sum(sim_img, axis = 0)
sim_sum = sum(sim_profile)

data_profile = np.sum(data_img, axis = 0)
x = np.linspace(0, len(data_profile), len(data_profile))

bp_index = posicion_segundo_maximo(data_profile)

plt.figure()
plt.plot(x, data_profile/max(data_profile), label = 'data' )
plt.plot(x + y_offcenter_px, sim_profile/max(sim_profile), label = 'sim')
plt.legend()





plt.figure()
plt.plot(x, data_profile, label = 'data' )
plt.legend()

plt.figure()
plt.plot(x + y_offcenter_px, sim_profile, label = 'sim')
plt.legend()

end = time.time()
print('Time elapsed:\t', end-start)