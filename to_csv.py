# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 17:41:56 2024

@author: usuario
"""

import numpy as np
import csv

P = np.array([65.6,15.0,6.4,3.5,2.25,1.55,1.05,0.81,0.61,0.49,0.39,0.3,
                 0.25,0.2,0.16,0.12,0.095,0.075,0.063])/100

n = np.array((range(1,20)))
with open("e_cl_distr.csv", "a+", newline ='') as csvfile:
    wr = csv.writer(csvfile, dialect='excel', delimiter=',')
    for i in range(len(P)):
        wr.writerow([n[i],P[i]])