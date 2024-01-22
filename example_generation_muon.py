# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:16:43 2023

@author: diego
"""

from Alpha_track_simulator import *
import numpy as np
import matplotlib.patches as patches
import scipy.special as spc

"""
Example of muon's generation with our simulation framework

In the following we are going to evaluate the ionization produced in a given 
medium (Argon, in this case) by muons of a given Energy. This will result in a 
number of clusters and electrons/cluster.
"""
plt.close('all')
# =============================================================================
# We take the Energy from a Poisson distribution centered in ~5 Gev. We can be 
# more precisse in this. (Reference in the main class code)
# =============================================================================

E = np.random.poisson(lam = 5000)   #E = 5000 [MeV]
muonAr = muon_ionization(E = E)
dEdx, n_cl_cm = muonAr.ionization(gas = 'Argon', exp_data = 'Argon') 

muons = []                          #list where muon tracks will be stored
dimensions = [10,10,15]
muon = muon_generator(energy = E, dEdx = dEdx, geometry = dimensions, n_cl_cm = n_cl_cm)
muon.produce_muon(n = 8, store = muons)


plt.figure()

plt.title('Muon clean tracks')
plt.plot([muons[0].z0, muons[0].z], [muons[0].y0, muons[0].y], color='r')
plt.plot([muons[1].z0, muons[1].z], [muons[1].y0, muons[1].y], color='g')
plt.plot([muons[2].z0, muons[2].z], [muons[2].y0, muons[2].y], color='b')
plt.plot([muons[3].z0, muons[3].z], [muons[3].y0, muons[3].y], color='b')
plt.plot([muons[4].z0, muons[4].z], [muons[4].y0, muons[4].y], color='b')
plt.plot([muons[5].z0, muons[5].z], [muons[5].y0, muons[5].y], color='b')
plt.plot([muons[6].z0, muons[6].z], [muons[6].y0, muons[6].y], color='b')
plt.plot([muons[7].z0, muons[7].z], [muons[7].y0, muons[7].y], color='b')
#Chamber
plt.vlines(x = 0, ymin = 0, ymax=dimensions[1], color='black' )
plt.vlines(x = dimensions[2], ymin = 0, ymax=dimensions[1], color='black' )
plt.hlines(y = 0, xmin = 0, xmax = dimensions[2], color = 'black')
plt.hlines(y = dimensions[1], xmin = 0, xmax = dimensions[2], color = 'black')
#####
plt.xlabel('z')
plt.ylabel('y')

"""
for muon in muons:
    
    muon.fill()
    
"""
muon1, muon2, muon3, muon4, muon5, muon6, muon7,muon8 = muons
zs1, ys1 = muon1.fill()
zs2, ys2 = muon2.fill()
zs3, ys3 = muon3.fill()
zs4, ys4 = muon4.fill()
zs5, ys5 = muon5.fill()
zs6, ys6 = muon6.fill()
zs7, ys7 = muon7.fill()
zs8, ys8 = muon8.fill()

fig, ax = plt.subplots()

h2d = ax.hist2d(x = zs1, y = ys1, bins = 100)
fig.colorbar(h2d[3], ax=ax)
#plt.show()


hist, x_edges, y_edges = h2d[0], h2d[1], h2d[2]

# Invertir el eje X
inverted_hist = np.flipud(hist.T)
inverted_x_edges = np.flip(x_edges)

plt.figure()
plt.imshow(inverted_hist, extent=[inverted_x_edges[-1], inverted_x_edges[0], y_edges[0], y_edges[-1]], cmap='viridis')
plt.xlabel('z [cm]', fontsize = 13)
plt.ylabel('y [cm]', fontsize = 13)
#plt.colorbar()


plt.tight_layout()
plt.show()


n = np.linspace(1,125,125)
def P_n_interac(r, D, lambd):
    
    return (1 / spc.factorial(r)) * (D / lambd)**r * np.exp(-D / lambd)

n_pdf = P_n_interac(n,muon1.track_length ,3.9170601906139784e-05)
"""
plt.figure()
#plt.hist(muon1.n_e_cl, bins = 125, density = True)
plt.plot(n,n_pdf, color = 'r')
"""