# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 12:15:34 2024

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
from Alpha_track_simulator import *
from scipy.optimize import curve_fit


plt.close('all')

#INPUTS========================================================================
n_tracks = 1
Emuon = 4000                    #(MeV) // assuming CR as source
muons = []                      #list to store muon objects
dimension = ([500,500,500])           #[x,y,z] (cm) of a rectangular-wise chamber
sigma_diff = 0.015
sigma_PSF = 0
pxs = 100

###############################################################################
#Pixel size calculation
px_size = dimension[0] / pxs       #cm / px


###############################################################################
"""
First of all we have to generate some muons so we use the simulator.
"""
print('Running...')
#Muons generation
muon = muon_generator(energy = Emuon, geometry = dimension)                    #First generate the muon object
muon.produce_muon(n = n_tracks, store = muons, line= True)                                 #generate muon's obj stored in muons list
#Applying diffusion to each track
diff_handler=Diffusion_handler(sigma_diff= sigma_diff ,sigma_PSF=sigma_PSF)
diff_handler.diffuse(muons)
#Filling tracks with e-
for muon in muons: muon.fill()
#Noise object
noise=Noise()
 

bins_x = int(muons[0].x * (1 / px_size))
bins_y = 100#int(5 * (1 / px_size))#int(muons[0].y  * (1 / px_size))

#Plot the final tracks
image2d=Image_2D(track_list = muons, hist_args={"bins":[bins_x, bins_y]})
image2d.track_plot()
image2d.plot_hist(noise_object=noise)
image2d.plot_x()
###############################################################################
e_2d = image2d.Hist2D       #e- positions in the 2d array (detector's pixels (?) )

total_e_col = []            #total electrons in a column of pixels
probab_centroid = []
total_e_col.append(np.sum(e_2d, axis= 1)) #CORREGIR!!!!
probab_centroid1 = []
for i in range(bins_y):
    if total_e_col[0][i] == 0:
        probab_centroid.append(np.zeros(len(total_e_col[0])))
    else:
        probab_centroid.append(e_2d[i] / total_e_col[0][i])          #FIXME

probab_centroid = np.array(probab_centroid)
probab_centroid = np.flip(probab_centroid, axis = 0)                 #Rotate the y-axis (numpy do weird things with the 2d arrays)


#print('probs:', probab_centroid)
plt.figure()
plt.imshow(probab_centroid, cmap = 'viridis')

"""
FROM HERE ON THERE COULD BE ERRORS 
"""
ML_position = []
for i in range(bins_y):
    ML_position.append(np.argmax(probab_centroid[i]))        #Most Likely position

ML_position = np.array(ML_position)
X = np.linspace(0, bins_y, bins_y)

plt.figure()
plt.plot(ML_position, X[::-1], 'o')

#LINEAR FIT====================================================================
def linear_fit(x, a, b):
    return a * x + b

fit_opt, fit_cov = curve_fit(linear_fit, ML_position, X[::-1])

plt.plot(ML_position, linear_fit(ML_position, *fit_opt))


full_px_array = np.zeros((pxs, pxs))
"""
start_x = int(muons[0].x0 * px_size)
start_y = int(muons[0].y0 * px_size)
# Insertar la imagen en la matriz de ceros de 100x100
full_px_array[start_y:start_y + len(probab_centroid), start_x:start_x + len(probab_centroid[0])] = probab_centroid
"""
"""
muon1 = muons[0]
if muon1.phi <= np.pi / 2:
    start_x = int(muons[0].x0 * (1/px_size))
    start_y = 0                 #int(muons[0].y0 * (1/px_size)) - len(probab_centroid[0])
    # Insertar la imagen en la matriz de ceros de 100x100
    full_px_array[start_y:start_y + len(probab_centroid), start_x:start_x + len(probab_centroid[0])] = probab_centroid

else:
    start_x = int(muons[0].x0 * (1/px_size))
    start_y = pxs - int(muons[0].y0 * (1/px_size))
    # Insertar la imagen en la matriz de ceros de 100x100
    full_px_array[start_y:start_y + len(probab_centroid), start_x:start_x + len(probab_centroid[0])] = probab_centroid

"""

plt.figure()
plt.imshow(full_px_array, cmap='viridis')
#plt.plot(ML_position, (50 + max(linear_fit(ML_position, *fit_opt))) - linear_fit(ML_position, *fit_opt), color = 'r')
plt.colorbar()
plt.show()





