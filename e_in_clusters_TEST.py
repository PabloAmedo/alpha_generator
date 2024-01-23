# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 16:24:39 2023

@author: usuario
"""

import numpy as np


# IONIZATION CALCULATION

"""
Should I implement ALL the cluster calculation here??? Like NIST data and stuff 
|
|--> Maybe as another method (?)
    |
    |-->class Ionization:
        methods: calculate the clusters per cm (gas)
        

Now I am taken the data from Clusters_calculation.py and moving on...
"""

n_clusters_length = 57.3            #[clusters / cm] 
mfp = 27.8                          #[1/cm] mean free path
n_avg = 3.3                         #[e/cl] average number of electrons per cluster
L = 150                             #[cm] length of the track (len of BAT (?))

n_cl_avg = np.random.poisson(lam = L * n_clusters_length)

def distr_1n2(n):
    """
    exp_prob = np.array([65.6,15.0,6.4,3.5,2.25,1.55,1.05,0.81,0.61,0.49,0.39,0.3,
                     0.25,0.2,0.16,0.12,0.095,0.075,0.063]) / 100
    if n < 19:
        return exp_prob[n]
    if n >= 19:
        return 1/n**2"""
    return 1/n**2

n_e_cl=[]
while len(n_e_cl) < n_cl_avg:
    a = np.random.randint(1,125)
    b = np.random.rand()
    
    if b < distr_1n2(a):
        n_e_cl.append(a)
        
#Computing totals, avg, and verifying
e_total = sum(n_e_cl)
cl_total = n_cl_avg
e_cl_avg = np.mean(n_e_cl)

#Showing results
print('Number of clusters (total):\t\t', cl_total)
print('Number of electrons (total):\t', e_total)
print('\n')
print('Number of e- per cluster (sim):\t', e_cl_avg)
print('Number of e- per cluster (teo):\t', n_avg)

"""
FIXME: (or update-me)

* Are Poisson distributions ok??
* Where the n_e_max goes here?? (I mean the 125 obtained in TESTSTUFF.py) 
"""
