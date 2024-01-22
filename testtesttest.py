# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 12:27:10 2023

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt

exp_prob = np.array([65.6,15.0,6.4,3.5,2.25,1.55,1.05,0.81,0.61,0.49,0.39,0.3,
                 0.25,0.2,0.16,0.12,0.095,0.075,0.063]) / 100

n_prob = np.linspace(19, 125, 106)

longarray = np.concatenate((exp_prob, (21.6/n_prob**2)/100), axis = 0)

n_exp = np.linspace(1,19,19, dtype = int)

plt.figure()

plt.title('n electrons distribution', fontsize = 20)
plt.plot(n_exp,exp_prob, 'x', color='black', label = 'experimental data [FIS91]' )
plt.plot(n_prob, (21.6/n_prob**2)/100, '--', color = 'r', label = 'extrapolation')
plt.xlabel('n electrons', fontsize = 13)
plt.ylabel('P', fontsize = 13)
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()