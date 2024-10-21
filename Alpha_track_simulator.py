# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 10:49:54 2023

@author: Pablo
"""
import numpy as np
import os as os
import pandas as pd
import scipy as scipy
import scipy.optimize as opt
import struct as struct
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import glob as glob
import math as math
import scipy.ndimage.filters as filters
from scipy.interpolate import interp1d

#Avoid warning in console
import warnings
from scipy.optimize import OptimizeWarning
warnings.filterwarnings("ignore", category = OptimizeWarning)                  #this avoids the covariance calculation warning 

import sys
current_path = os.getcwd()

auxiliar_path = current_path + '\\auxiliar_scripts'
examples_path = current_path + '\\examples'
data_path = current_path + '\\data'
data_path2 = current_path + '\\data\\clusters_distributions'

add_to_path = [auxiliar_path,examples_path, data_path,  data_path2]

for path in add_to_path:
    if add_to_path in sys.path:
        break
    else:
        sys.path.append(path)

from dEdx_tools import *
from general_tools import *
from bragg_distribution import *



class Source:
    
    """This is the source class. It includes all the relevant information about the source"""
    
    def __init__(self, rate = 500, energy = 5.5, radius = 0.35, M = 933, range_alpha = 4.8, We = 25.7e-6, red_fact = 1):
        
        self.rate=rate #In Hz
        
        self.energy=5.5 #In Mv
        
        self.radius=radius # In cm
        
        self.range_alpha=range_alpha
        
        self.M=M
        self.z=2
        self.n_e=int(energy/We)                                  #Add We at some point
        self.red_fact = red_fact
    
    def produce_alpha(self, n, store,ionization_profile, phi_in=None, ath_in=None, theta=None): 
        """This method produces an alpha track from the alpha_tracks class
        
            - n                 : number of alpha tracks
            - store             : list to store results
            - ionization_profile: mainly Bragg distribution for ionization
            - phi_in            : angle in xy plane [0, 2*pi]
            - ath_in            : angle in yz plane [0, pi/2]
            - (!)theta          : angle inside source active region --> determine initial position (x0, y0)
        
        """        
        n_electrons = int(214007/self.red_fact)
        if ionization_profile=='Bragg':
            self.x,self.Sp,self.acum=bragg_peak(self.energy,self.range_alpha)
        else:
            self.x=None
            self.Sp=None
            self.acum=None

        for i in range(n):
            if phi_in == None:
                phi = np.random.rand() * 2 * np.pi
            else:
                phi = phi_in
                
            if ath_in == None:
                ath = np.random.rand() * np.pi / 2
            else:
                ath = ath_in
                            
            theta = np.random.rand() * 2 * np.pi

            x0 = np.cos(theta) * self.radius * np.random.rand()
            y0 = np.sin(theta) * self.radius * np.random.rand()
                 
            alpha = Alphas_tracks(x = self.x, acum = self.acum, range_alpha = self.range_alpha,
                                phi = phi, ath = ath, x0 = x0, y0 = y0, ionization_profile = ionization_profile, 
                                n_electrons = n_electrons)
            store.append(alpha)
            
        return store
    
    def collimator(self, stored_tracks, collimated_tracks, phi_col, slit, dist):
        
        """
        Method to set a collimation for the alphas. 
        
        stored_tracks       --> list with tracks previously generated
        collimated_tracks   --> list where the collimated tracks will be stored
        phi                 --> collimation angle
        slit                --> width of the window
        dist                --> distance from the alpha to the window
        """
        # beta. Angle between source and window's edge
        beta = np.arctan((slit / 2 / dist))
        
        for track in stored_tracks:
            if track.phi <= phi_col + beta and track.phi >= phi_col-beta:
                collimated_tracks.append(track)
        
        return collimated_tracks



class muon_generator:#FIXME: change name to cp_generator -- charged particle
    
    """
    muon_generator class will generate 'n' muon tracks 
    with the corresponding number of ionization clusters + e-. 
    -----------------------------------------------------------------------
        -energy            : (float)                energy in MeV
        -geometry          : (list/array of float)  [x,y,z] dimensions of the chamber
        -mass              : (float)                mass in MeV/c2
        -gas               : (str)                  options: 'Argon', 'CH4', 'ArCH4-90/10'
        -pressure          : (float)                pressure in bar
    -----------------------------------------------------------------------
    Tracks are generated randomly in a choosen rectanglar volume. We make the 
    next assumptions:
        -We are considering that muon flux is coming from one lateral only.
        -We are only considering muons that go throught the entire region.
        -Not considering muons generated inside due to interactions.
        -Not considering muons stopping in the chamber region.
        -Muons do not loss energy
        
    """
    
    def __init__(self, energy = 5e3,  geometry = [1,1,1], mass = 105.66, gas = 'Argon', pressure = 1, red_fact = 1):
        #Define the initial conditions here
        self.energy = energy                                                   #MeV
        self.xmax, self.ymax, self.zmax = geometry
        self.mass = mass                                                       #MeV/c2 muon mass by default
        self.gas = gas
        self.P = pressure                                                      #bar
        
        
    def produce_muon(self, n, store, position_in = None, phi_in = None, ath_in = None, line = False, 
                     e_cut = 10000, n_cl_cm_in = None, length = None):
        
        """
        This method is used to generate n muon tracks by randomly generate a 
        y0 position where the muon will enter the chamber and also phi0 and 
        theta0 defining the direction the muon will follow across it.

        -----------------------------------------------------------------------
            -n                 : (int) number of tracks to generate
            -store             : (list / array-like) to store tracks
            -y0                : (float) 
            -phi0              : (float)  initial position
            -theta0            : (float) 
            -e_cut             : (int) max number of e- produced in a single cl
            -(!)line              : (bool) trackas crossing in a straight line the chamber
        -----------------------------------------------------------------------
        
        Output: muon_tracks object
        """
        
        ########################### POSITIONS #################################
        for i in range(n):
            if position_in:
                self.x0, self.y0, self.z0 = position_in
            else:
                #If initial position is not given we generate a random one
                self.x0 = 0 # np.random.rand(n) * self.xmax     # IF X0 = 0 THE TRACK 'COMES' FROM THIS SIDE
                self.y0 = np.random.rand() * self.ymax
                self.z0 = np.random.rand() * self.zmax
            
            #Get/generate angles
            self.phi = phi_in if phi_in != None else np.random.rand() * np.random.choice((-1, 1)) * np.pi / 2 
            self.ath = ath_in if ath_in != None else np.random.rand() *  np.pi
            
         
            phi_min = np.arctan((self.ymax - self.y0 ) / self.xmax)
            ath_min = np.arctan(self.xmax / (self.zmax - self.z0 ))
            ath_max = np.arctan(self.z0 / self.xmax) + np.pi / 2
    
            #Computing final positions
            if length:
                xout = length * np.cos(self.phi) * np.cos(self.ath)
                yout = length * np.cos(self.phi) * np.sin(self.ath)
                zout = length * np.cos(self.phi) #FIXME 
                    
            else:
                #If a length is not given we assume it goes through the entire region
                xout = self.xmax if self.phi <= phi_min else (self.ymax - self.y0) / np.tan(self.phi)
                yout = self.y0 + np.tan(self.phi) * xout if self.phi <= phi_min else self.ymax 
                
                if self.ath <= ath_min:
                    zout = self.zmax
                elif ath_min < self.ath <= ath_max:
                    zout = self.z0 + abs(xout - self.x0) * np.tan(np.pi/2 - self.ath) if self.ath <= np.pi/2 else self.z0 + abs(xout - self.x0) * np.tan(self.ath - np.pi/2)
                else:
                    zout = 0
            tr_len = np.sqrt((xout - self.x0)**2 + (yout - self.y0)**2 + (zout - self.z0)**2)
        
            ####################### CLUSTER + ELECTRONS #######################
            #Getting the number of clusters for the given energy and particle (self.mass)
            n_cl_cm = self.P * dNdx(self.energy, self.mass, pure_argon = False) #if n_cl_cm_in == None else self.P * n_cl_cm_in #from HEED simulations
            n_cl_avg = np.random.poisson(lam = tr_len * n_cl_cm)
            
            #FIXME: I HAVE TO CLEAN THIS A BIT (remove the useless ones)
            #CLUSTER DISTRIBUTOIN FOR DIFFERENT GASES --- Maybe this could be done from a dict???? - fixme
            if self.gas == 'Argon':#Fischle 
                P_el = np.loadtxt('data/clusters_distributions/Ar_cluster_distribution_experimental.txt') / 100
                n_el = np.linspace(1, 19, 19)
                probabilities_20 = ClusterParametrizationAr(np.linspace(20, e_cut, e_cut - 19))
                probs_final = np.concatenate((P_el, probabilities_20))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
            elif self.gas == 'CH4':#Fischle
                P_el = np.loadtxt('data/clusters_distributions/CH4_cluster_distribution_experimental.txt') / 100
                n_el = np.linspace(1, 19, 19)
                probabilities_20 = ClusterParametrizationCH4(np.linspace(20, e_cut, e_cut - 19))
                probs_final = np.concatenate((P_el, probabilities_20))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                
            elif self.gas == 'ArCH4_93-7':#from HEED 
                data_ArCH4_9307 = np.loadtxt('data/clusters_distributions/ArCH4-93-7-10bar_0.1x0.1x0.1mmCell_pi2.5GeV_normalized.txt' , skiprows = 1)
                n_el, P_el, _ = np.split(data_ArCH4_9307, 3, axis = 1)
                P_el = P_el.flatten() ; n_el = n_el.flatten()
                probabilities_aboveAr = ClusterParametrizationAr(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probabilities_aboveCH4 = ClusterParametrizationCH4(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probabilities_above = probabilities_aboveAr*0.93 + probabilities_aboveCH4*0.07
                probs_final = np.concatenate((P_el, probabilities_above))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                
            elif self.gas == 'ArCH4_90-10':#from HEED 
                data_ArCH4_9010 = np.loadtxt('data/clusters_distributions/ArCH4-90-10-10bar_0.1x0.1x0.1mmCell_pi2.5GeV_normalized.txt' , skiprows = 1)
                n_el, P_el, _ = np.split(data_ArCH4_9010, 3, axis = 1)
                P_el = P_el.flatten() ; n_el = n_el.flatten()
                probabilities_aboveAr = ClusterParametrizationAr(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probabilities_aboveCH4 = ClusterParametrizationCH4(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probabilities_above = probabilities_aboveAr*0.90 + probabilities_aboveCH4*0.1
                probs_final = np.concatenate((P_el, probabilities_above))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                
                
            elif self.gas == 'ArCF4_99-1':#from HEED 
                data_ArCF4_9901 = np.loadtxt('alpha_generator/data/clusters_distributions/ArCF4-99-1-10bar_0.1x0.1x0.1mmCell_8GeVmuons_normalized.txt' , skiprows = 1)
                n_el, P_el, _ = np.split(data_ArCF4_9901, 3, axis = 1)
                P_el = P_el.flatten() ; n_el = n_el.flatten()
                probabilities_above = ClusterParametrizationAr(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probs_final = np.concatenate((P_el, probabilities_above))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                
                
            elif self.gas == 'DEGRAD_data': #DEGRAD (it was used to compare with Fischle data)
                data_DEGRAD = np.loadtxt('data/clusters_distributions/Degrad_Pure_Ar_v3.txt')
                n_el, P_el = np.split(data_DEGRAD, 2, axis = 1)
                P_el = P_el.flatten() ; n_el = n_el.flatten()
                probabilities_above = ClusterParametrizationAr(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probs_final = np.concatenate((P_el, probabilities_above))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())

                
            elif self.gas == 'PEP4':#from HEED (01mm to remove the XR propagation)
                n_cl_cm = 0.9844562192883783 * self.P * dNdx(self.energy, self.mass)
                n_cl_avg = np.random.poisson(lam = tr_len * n_cl_cm)
                
                data_PEP4 = np.loadtxt('data/clusters_distributions/ArCH4-80-20-10bar_0.1x0.1x0.1mmCell_2.5GeVmuons_normalized.txt' , skiprows = 1)
                n_el, P_el, _ = np.split(data_PEP4, 3, axis = 1)
                P_el = P_el.flatten() ; n_el = n_el.flatten()
                probabilities_aboveAr = ClusterParametrizationAr(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probabilities_aboveCH4 = ClusterParametrizationCH4(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probabilities_above = probabilities_aboveAr*0.80 + probabilities_aboveCH4*0.20
                probs_final = np.concatenate((P_el, probabilities_above))
                
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                
            elif self.gas == 'nano-N2':
                #n_cl_cm = 2 #DELETE
                #n_cl_avg = np.random.poisson(lam = tr_len * n_cl_cm)
                data_nanoN2 = np.loadtxt('nanodosimetry/cluster_size/nitrogen-10bar_0.1x0.1x0.1mmCell_2.5GeVmuons_normalized.txt')
                n_el, P_el, _ = np.split(data_nanoN2, 3, axis = 1)
                P_el = P_el.flatten() / 100 ; n_el = n_el.flatten()
                #probabilities_above = ClusterParametrizationGeneral(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probs_final =P_el[:e_cut] # np.concatenate((P_el, probabilities_above))
                
                print('a=', n_el[:e_cut])
                print('p=', len(probs_final / probs_final.sum()))
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                
            elif self.gas == 'nano-C3H8':
                
                n_cl_avg = np.random.poisson(lam = tr_len[i] * n_cl_cm * 2.65)
                data_nanoC3H8 = np.loadtxt('nanodosimetry/cluster_size/propane-10bar_0.1x0.1x0.1mmCell_2.5GeVmuons_normalized.txt')
                n_el, P_el, _ = np.split(data_nanoC3H8, 3, axis = 1)
                P_el = P_el.flatten() / 100 ; n_el = n_el.flatten()
                #probabilities_above = ClusterParametrizationGeneral(np.linspace(int(max(n_el))+1, e_cut, e_cut-int(max(n_el))))
                probs_final = P_el[:e_cut] # np.concatenate((P_el, probabilities_above))
                
                print('a=',e_cut)
                print('p=', len(probs_final))# / probs_final.sum()))
                n_e_cl = np.random.choice(np.linspace(1, e_cut, e_cut), size = n_cl_avg, 
                                          p = probs_final / probs_final.sum())
                

            #Uniform distribution of the clusters along the track
            cl_position = np.random.uniform(low = 0, high = xout, size = len(n_e_cl)) #x position
            cl_position = np.sort(cl_position)
            
            #Saving into muon_tracks object
            muon = muon_tracks(energy = self.energy, track_length = tr_len, initial_position = [self.x0, self.y0, self.z0], 
                               final_position =[xout, yout, zout] , ath = self.ath, phi= self.phi,
                               clusters = cl_position, n_e_cl = n_e_cl, geometry = [self.xmax, self.ymax, self.zmax])
            store.append(muon)
            
        return n_e_cl, cl_position


class muon_tracks(muon_generator):
    
    """
    muon_tracks class. It contains all the relevant information about the 
    charged particles traveling across the chamber.
    
    *This has been copied form alpha_tracks almost entirely! To be revised
    """
    
    def __init__(self, energy, track_length, initial_position, final_position, ath, phi,  clusters, n_e_cl, geometry = [1,1,1], spread = 0):
        
        self.energy = energy
        self.track_length = track_length
        self.x0, self.y0, self.z0 = initial_position
        self.x, self.y, self.z = final_position
        self.ath = ath
        self.phi = phi
        self.clusters = clusters    #x coord of each cluster
        self.n_e_cl = n_e_cl
        self.n_electrons = np.sum(n_e_cl)
        self.xmax, self.ymax, self.zmax = geometry
        
        if self.y > self.y0:
            self.theta = np.arctan(abs(self.y-self.y0) / (self.x-self.x0))
        elif self.y == self.y0:
            self.theta = np.pi / 2
        else:
            self.theta = np.pi / 2 + np.arctan(abs(self.y-self.y0) / (self.x-self.x0))
        
        #This is defined at 2 sigma and will be handled by the Difussion_handler class
        self.spread=spread

    
    
    def fill(self, diff = False):
        """
        This method is used to fill the clean muon track through the chamber 
        with the electrons produced in the ionization.
        """
        #We have the initial and final coordinates in XY. Assuming a linear track
        #we have to calculate the slope to have the linear equation
        #y(x) = m * x + n
        
        slope_xy = (self.y-self.y0) / (self.x-self.x0)
        slope_xz = (self.z-self.z0) / (self.x-self.x0)
        linear_track_xy = lambda x : slope_xy * x + self.y0
        linear_track_xz = lambda x : slope_xz * x + self.z0
        cl_y_tr = []                  #List to store y coord of cl around track 
        cl_x_tr = []                  #List to store x coord of cl around track
        cl_z_tr = []
        ne = []
        
        for i in range(len(self.clusters)):
            cl_y_pos = linear_track_xy(self.clusters[i])
            cl_z_pos = linear_track_xz(self.clusters[i])
            if (0 <= cl_y_pos <= self.ymax) and (0 <= cl_z_pos <= self.zmax):
                cl_y_tr.append(float(cl_y_pos))
                cl_x_tr.append(self.clusters[i])
                cl_z_tr.append(float(cl_z_pos))
                ne.append(self.n_e_cl[i])
        
        cl_x_tr = np.array(cl_x_tr)
        cl_y_tr = np.array(cl_y_tr)
        cl_z_tr = np.array(cl_z_tr)
        
        #Storage of electrons (z,y) positions
        self.electron_positions=np.zeros([int(sum(ne)), 3]) # coordenadas [x,y] para los 50 electrones
        #Storage of electrons (z,y) positions after diffusion
        self.electron_positions_diff=np.zeros([int(sum(ne)), 3])
        
        #Assign positions to e in cluster 
        aux_electrons_total = 0
        for i in range(len(cl_y_tr)):#FIXME
            for j in range(int(self.n_e_cl[i])):
                #Change the positions
                self.electron_positions[aux_electrons_total, 0] = cl_x_tr[i]
                self.electron_positions[aux_electrons_total, 1] = cl_y_tr[i]
                self.electron_positions[aux_electrons_total, 2] = cl_z_tr[i]
                
                #Change the positions
                self.electron_positions_diff[aux_electrons_total, 0] = cl_x_tr[i]
                self.electron_positions_diff[aux_electrons_total, 1] = cl_y_tr[i]
                self.electron_positions_diff[aux_electrons_total, 2] = cl_z_tr[i]
                
                aux_electrons_total += 1

        # Calculate the transverse diffusion on the z-y plane
        if diff == True:
            #Now we need to draw the number from the gaussian distribution
            pos = np.random.multivariate_normal((0,0,0), cov = self.spread**2 * np.identity(3), size = int(sum(ne)))
            
            self.electron_positions_diff[:,0] += pos[:,0]
            self.electron_positions_diff[:,1] += pos[:,1]
            self.electron_positions_diff[:,2] += pos[:,2]
            
            #return cl_x_tr, cl_y_tr
    
    
    def select(self,variable_scan="x",min_var=0,max_var=100,array_to_store=None):#FIXME --> Change for loop limits!
        """This method scans the positions of the electrons either in the x or y direction and
        it selects them based on a condition on the other variable. Ex: we scan in x and select
        electrons whose y position fulfills y_min<y<y_max"""
        
        #Select the apropiate index, the opossite of variable_scan
        if variable_scan=="x":
            index=1;other_index=0
        else:
            index=0;other_index=1
            
        #Iterate over the electrons
        for j in range(int(self.n_electrons)):
            #Check if they are whitin our current 
            if min_var<self.electron_positions_diff[j,index]<=max_var:
                array_to_store.append(self.electron_positions_diff[j,other_index])
                
        return array_to_store
    
        
class Alphas_tracks(Source):
    
    """This is the alpha class. It contains all the relevant information about the alpha track,
    like the ionization profile, positions of electrons, etc"""
    
    def __init__(self,x,acum,range_alpha,phi=None,ath=None,x0=None,y0=None,spread=0,ionization_profile="Flat", n_electrons=int(5.5/25.7e-6)): #FIXME
        super().__init__()
        self.ionization = ionization_profile
        self.X = x
        self.acum  = acum            #si no da error
                                  #son vectores x y Sp, necesarios para los randoms
        
        #This is the range of the alpha. Get it from NIST
        self.range_max = range_alpha #In cm
        
        #Athimutal angle. This reduces our effective alpha range by the projection
        self.ath = ath
        
        #Define an effective range
        self.range = self.range_max * np.cos(ath)
                      
        self.x0 = x0
        self.y0 = y0
      
        #X position in the plane
        self.x = np.cos(phi) * self.range * np.cos(ath) + self.x0 
        
        #Y position in the plane
        self.y = np.sin(phi) * self.range * np.cos(ath) + self.y0
        
        #Phi angle
        self.phi = phi
               
        #This is defined at 2 sigma and will be handled by the Difussion_handler class
        self.spread = spread
        
        #Set the ionization profile to be flat between the 0 and the range
        self.dict_ion_prof = {
            "Flat": np.random.rand,
            "Bragg": random_bragg
            }
        
        self.ionization_profile=self.dict_ion_prof[self.ionization]
        
        #Number of electrons in a track
        # self.n_electrons=source.energy/Gas.I
        # self.n_electrons=50
        self.n_electrons = n_electrons / self.red_fact
        
        #Storage of electrons (x,y) positions
        self.electron_positions=np.zeros([n_electrons,2]) # coordenadas [x,y] para los 50 electrones
        
        #Storage of electrons (x,y) positions after diffusion
        self.electron_positions_diff=np.zeros([n_electrons,2])
        
    def fill(self, diff = False):
        
        """This method must only be called when we want to fill the final alpha track with
        electrons. 
        
        First it picks a random number following a given ionization profile to situate the electron
        along the track. Then it picks a number to situate along transverse gaussian distribution
        that mimics diffusion. If the spread is zero then we don't apply the transversal diffusion.
        """

        #Draw a number from the ionization profile distribution
        if np.all(self.electron_positions == 0):
            for i in range(int(self.n_electrons)):
                if self.ionization == 'Bragg':
                    r_pos = self.ionization_profile(self.X, self.acum) * self.range
                else:
                    r_pos = self.ionization_profile() * self.range
                
                self.electron_positions[i,0] = np.cos(self.phi) * r_pos + self.x0
                self.electron_positions[i,1] = np.sin(self.phi) * r_pos + self.y0
                
                self.electron_positions_diff[i,0] =  np.cos(self.phi) * r_pos + self.x0
                self.electron_positions_diff[i,1] =  np.sin(self.phi) * r_pos + self.y0

        # Calculate the transverse diffusion on the x-y plane
        if diff == True:
            #Now we need to draw the number from the gaussian distribution
            pos = np.random.multivariate_normal((0,0), cov = self.spread**2 * np.identity(2), size = int(self.n_electrons))
            
            self.electron_positions_diff[:,0] += pos[:,0]
            
            self.electron_positions_diff[:,1] += pos[:,1]
                
                
    def select(self, variable_scan = "x", min_var = 0, max_var = 100, array_to_store = None):
        """This method scans the positions of the electrons either in the x or y direction and
        it selects them based on a condition on the other variable. Ex: we scan in x and select
        electrons whose y position fulfills y_min<y<y_max"""
        
        #Select the apropiate index, the opossite of variable_scan
        if variable_scan=="x":
            index=1;other_index=0
        else:
            index=0;other_index=1
            
        #Iterate over the electrons
        for j in range(int(self.n_electrons)):
            #Check if they are whitin our current 
            if min_var<self.electron_positions_diff[j,index]<=max_var:
                array_to_store.append(self.electron_positions_diff[j,other_index])
                
        return array_to_store


class Gas:
    
    """This is the gass class. It includes information about the W, drift velocity, etc"""
    
    def __init__(self,gas = 'ArCF4_99-1', vz = None, Dt = None, Dl = None,
                 Wi = 25.7,density=1.662e-3, A = 40, Z = 18, I = 188e-6,
                 L_drift = 14.2, Pressure = 1):
        
        self.gas = gas
        
        #This has to be changed for the case of a mixture
        #Default values are for Pure Ar
        self.Wi = Wi                #eV/e
        self.density = density      #g/cm3
        self.A = A
        self.Z = Z      
        self.I = I                  #Mev
        
        #Setup info --> change this
        self.L_drift = L_drift
        self.Pressure = Pressure
        
        
    def vz(self, RedField):
        #Load data from: https://github.com/UTA-REST/ArXe_plus_anything/blob/master/Argon_Plots/Argon_CF4.png
        data_vz = pd.read_csv('alpha_generator/data/diffusion/driftV_' + self.gas + '.csv', delimiter = ';')
        #Linear interpolator to get vz for every Drift Field value
        vz_interpolate = interp1d(data_vz['Reduced Drift Field (V/cm/bar)'], data_vz['Drift Velocity (mm/us)'])
        self.drift_velocity = vz_interpolate(RedField)                          #value in mm/us
        
    def Dt(self, RedField):
        
        #Load data from: https://github.com/UTA-REST/ArXe_plus_anything/blob/master/Argon_Plots/Argon_CF4.png
        data_Dt = pd.read_csv('alpha_generator/data/diffusion/diffT_' + self.gas + '.csv', delimiter = ';')
        #Linear interpolator to get Dt for every Drift Field value
        Dt_interpolate = interp1d(data_Dt['Reduced Drift Field (V/cm/bar)'], data_Dt['Dt (sqrt(bar)*um*1/sqrt(cm))'])
        self.Dt_coef = Dt_interpolate(RedField)
        #Sigma diff (transversal) in mm
        print('Dtcoef', self.Dt_coef)
        self.sigma_diff = self.Dt_coef / np.sqrt(self.Pressure) * np.sqrt(self.L_drift) * 1e-3    #value in mm
        
    def Dl(self, RedField):
        """Currently this is not being used"""
        #Load data from: https://github.com/UTA-REST/ArXe_plus_anything/blob/master/Argon_Plots/Argon_CF4.png
        data_Dl = pd.read_csv('alpha_generator/data/diffusion/diffL_' + self.gas + '.csv', delimiter = ';')
        #Linear interpolator to get Dl for every Drift Field value
        Dl_interpolate = interp1d(data_Dl['Reduced Drift Field (V/cm/bar)'], data_Dl['Dl (sqrt(bar)*um*1/sqrt(cm))'])
        self.Dl_coef = Dl_interpolate(RedField)
        #Sigma diff (transversal) in mm
        self.sigma_diffL = self.Dl_coef / np.sqrt(self.Pressure) * np.sqrt(self.L_drift) * 1e-3    #value in mm

class Diffusion_handler(Gas):

    """This class handles the application of diffusion from the gas. It inherits from the gas
    class."""
    
    def __init__(self, sigma_diff = None, sigma_PSF = 0):
        
        super().__init__()
        if sigma_diff:
            self.sigma_diff = sigma_diff
        
        self.u = 1
        self.sigma_PSF = sigma_PSF
        self.sigmas = []
    def diffuse(self,alpha_list):
        
        #Take an alpha track and add a difussion in cm
        for track in alpha_list:
            track.spread= ( self.sigma_diff**2 + self.sigma_PSF )**0.5  # Desviación estándar de la Gaussiana (7.5 mm??? too much???)
    
    def diffuseL(self, track_list, RedField):

        #Take an alpha track and add a difussion in cm
        for track in track_list:
            #FIXME To get electron position you have to use fill method
            track.fill(diff = False)
            y_positions = track.electron_positions[:,1]
            self.z_positions = abs(y_positions) * np.tan(track.ath)

            ldrift = self.L_drift - self.z_positions
            
            for i in range(len(ldrift)):
                self.L_drift = ldrift[i]
                self.Dt(RedField)
                #print(ldrift[i])
                self.sigmas.append(self.sigma_diff)
            track.spread = ( (self.sigma_diff / 10)**2 + self.sigma_PSF**2 )**0.5
    
    
class Noise():

    """This class handles the noise addition to the final 2D histogram where we assume we have
    our array of pixels
    
    We care about the spread of the noise, not so much about the absolute value, so this noise
    magnitudes will probably be about spread
    
    """
    
    def __init__(self,dark_current=109):
        
        #Dark_current (electrons/s) of the camera is 190 e/s. This depends exponentially on the
        #temperature
        self.dark_current = dark_current
    
    def add_noise(self, exposition_time, Hist2d, sigma=5.7):
        
        """This method adds noise to each pixel of a 2D histogram as a function of the
        electronic noise and the exposure time"""
        
        #Electronic noise (e) = dark_current (e/s) * exposition_time (s)
        #These units are in electrons
        electronic_noise = self.dark_current #* exposition_time
        
        #Create an array of n_bins*n_bins and populate them with random numbers
        #random_gauss=np.random.normal(scale=electronic_noise,size=Hist2d.shape)
        random_gauss=np.random.normal(loc=electronic_noise,size=Hist2d.shape,scale=sigma)
        #Ignore the negative ones?
        #random_gauss=abs(random_gauss) #!!!!
        
        #Update the 2d histogram with this values and return it
        
        return Hist2d+random_gauss
    
    def add_noise_abs(self,Hist2d,mean=0,sigma=5.7):
        
        """This method adds noise to each pixel of a 2D histogram using a gaussian distribution
        for each pixel"""
        
        random_gauss=np.random.normal(loc=mean,size=Hist2d.shape,scale=sigma)
        #Ignore the negative ones?
        random_gauss=abs(random_gauss)
        
        #Update the 2d histogram with this values and return it
        
        return Hist2d + random_gauss
        
        

class Image_2D():
    
    """This class is essentially a 2D histogram that also stores extra information at the truth level:
        how many tracks were produced, their positions, the positions of the electrons, etc. This
        is done with numpy.histogram2d and the axis are already inverted..
        
        This object will be exported and loaded by our analysis framework"""
        
    def __init__(self,track_list, hist_args = {"bins":100}, QE = 1, GE = 1, Tp = 1, gain = 1, ph_poiss = False, pixel_size = 0.2):
        
        #List of tracks
        self.track_list = track_list
        
        #Get a list of all electron positions in all tracks
        #self.x_pos = np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,0] for i in range(len(track_list))]))
        #self.y_pos = np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,1] for i in range(len(track_list))]))
        self.x_pos = np.concatenate([track_list[i].electron_positions_diff[:,0] for i in range(len(track_list))])
        self.y_pos = np.concatenate([track_list[i].electron_positions_diff[:,1] for i in range(len(track_list))])
        
        #Get a list of all the true electron's positions in all tracks
        #self.x_pos_true = np.ndarray.flatten(np.asarray([track_list[i].electron_positions[:,0] for i in range(len(track_list))]))
        #self.y_pos_true = np.ndarray.flatten(np.asarray([track_list[i].electron_positions[:,1] for i in range(len(track_list))]))
        
        self.x_pos_true = np.concatenate([track_list[i].electron_positions[:,0] for i in range(len(track_list))])
        self.y_pos_true = np.concatenate([track_list[i].electron_positions[:,1] for i in range(len(track_list))])
        
        
        
        #Define the quantum efficiency, the geometrical efficiency, the transparency and the gain
        self.QE = QE
        self.GE = GE
        self.Tp = Tp
        self.gain = gain    
        """
        #hist_args['bins'] are the TOTAL number of pixels per axis for the entire detector.
        #We have to select the number of pixels involved in the tracking for smaller tracks.
        tlen_x = []
        tlen_y = []
        for track in track_list:
            tlen_x.append(track.x)
            tlen_y.append(abs(track.y - track.y0))
        
        track_length_x = max(tlen_x)
        track_length_y = max(tlen_y) if max(tlen_y) != 0.0 else 1.4
        binsx = int(track_length_x / pixel_size)
        binsy = int(track_length_y / pixel_size)
        
        hist_args['bins'] = (abs(int(binsx)), abs(int(binsy)))
        #Get the 2d hist for all the electrons and the edges from the list of tracks.
        #Extra arguments can be passed to the 2D numpy histogram function
        self.Hist2D_e,self.x_edges,self.y_edges=np.histogram2d(x = self.x_pos, y = self.y_pos, bins = (abs(int(binsx)), abs(int(binsy))))
        #Since numpy inverts the (y,x) histogram, invert it again
        self.Hist2D_e=self.Hist2D_e.T"""
        
        #Get the 2d hist for all the electrons and the edges from the list of tracks.
        #Extra arguments can be passed to the 2D numpy histogram function
        self.Hist2D_e, self.x_edges, self.y_edges = np.histogram2d(x = self.x_pos, y = self.y_pos, **hist_args)
        #Since numpy inverts the (y,x) histogram, invert it again
        self.Hist2D_e = self.Hist2D_e.T
        
        if ph_poiss==True:
            #We generate in each pixel a poisson distribution with lamda=nº of electrons in that
            #pixel
            self.Hist2D =np.random.poisson(self.Hist2D_e*self.gain)*self.QE*self.GE*self.Tp
            
        else:
            #Apply the quantum efficiency, the geometrical efficiency, the transparency and the gain to th
            #histogram. This is the fast way of getting the number of photons for each pixel.
            self.Hist2D = self.Hist2D_e*self.QE*self.GE*self.Tp*self.gain
    
    def track_plot(self,fig_in=None,axis_list_in=None):
        #This does the first plot (true tracks + electrons positions)
        """This plots the generated tracks in a 2D image. It's advisable not to do it if the number
        of tracks is too high because it will take time"""
        
        #If a figure and axis are provided, use them. Otherwise come up with our own
        if fig_in==None or axis_list_in==None:
            fig=plt.figure(figsize=(6,6))
            ax=fig.add_subplot(1,2,1)
            ax2=fig.add_subplot(1,2,2)
            
        else:
            fig=fig_in
            ax=axis_list_in[0]
            ax2=axis_list_in[1]
        
        #===Plot section===
        #This is super slow if you have to do it track by track. Just skip it if there are too many
        if len(self.track_list)<100:
            
            #Plot the original tracks
            for i in range(len(self.track_list)):
                
                ax.plot([self.track_list[i].x0,self.track_list[i].x],[self.track_list[i].y0,self.track_list[i].y])        
        
            ax.set_xlabel("x (cm) ")
            ax.set_ylabel("y (cm) ")
            ax.set_title("Original tracks")
        
            #Plot the tracks
            for i in range(len(self.track_list)):
                
                ax2.scatter(self.track_list[i].electron_positions[:,0],self.track_list[i].electron_positions[:,1],marker="o", color='r', zorder=3)
        
            #Restart the color cyle of the axis
            ax2.set_prop_cycle(None)
            
            #Plot the tracks with the diffused electrons
            for i in range(len(self.track_list)):
                
                ax2.scatter(self.track_list[i].electron_positions_diff[:,0],self.track_list[i].electron_positions_diff[:,1],marker="x")
            
            ax2.set_xlabel("x (cm) ")
            ax2.set_ylabel("y (cm) ")
            ax2.set_title("Electron positions")
            
        return fig
            
    def plot_hist(self, fig_in = None, axis_list_in = None, noise_object = None, exposition_time = 0):
        """This simply plots the 2D histogram"""
    
        H, yedges, xedges = self.Hist2D, self.y_edges, self.x_edges
        
        #If a figure and axis are provided, use them. Otherwise come up with our own
        if fig_in == None or axis_list_in == None:
            #Create a figure
            fig, (ax1,ax2) = plt.subplots(ncols=2)
            
        else:
            fig = fig_in
            ax1 = axis_list_in[0]
            ax2 = axis_list_in[1]
            
            
        #Plot it
        ax1.pcolormesh(xedges, yedges, H, cmap = 'rainbow',vmin = np.min(H), vmax = np.max(H))
        
        #ax1.plot(x, 2*np.log(x), 'k-')
        
        ax1.set_xlim(self.x_pos.min(), self.x_pos.max())
        ax1.set_ylim(self.y_pos.min(), self.y_pos.max())
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_title('Nº of tracks '+ str(len(self.track_list))+" , equivalent to "+str(len(self.track_list)/500)+" s exposure time")
        
        #Add the noise
        H2=noise_object.add_noise(exposition_time=exposition_time, Hist2d=H)
        
        #Plot it
        mesh2=plt.pcolormesh(xedges, yedges, H2, cmap='rainbow',vmin = np.min(H), vmax = np.max(H))
        ax2.pcolormesh(xedges, yedges, H2, cmap='rainbow',vmin = np.min(H), vmax = np.max(H))
        
        fig.colorbar(mesh2, ax=[ax1, ax2])
        
        #ax1.plot(x, 2*np.log(x), 'k-')
        
        ax2.set_xlim(self.x_pos.min(), self.x_pos.max())
        ax2.set_ylim(self.y_pos.min(), self.y_pos.max())
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_title('Nº of tracks '+ str(len(self.track_list))+" , equivalent to "+str(len(self.track_list)/500)+" s exposure time with noise")
        
        return fig
    
    def plot_x(self, variable_cut = 'x', min_var = -0.5, max_var = 0.5):
        electron_cut = []
        
        if variable_cut == 'x':
            other_variable = 'y'
        else:
            other_variable = 'x'
            
        for track in self.track_list:
            electron_cut = track.select(variable_cut, min_var = min_var, max_var = max_var, array_to_store = electron_cut)
        
        le_hist = np.histogram(electron_cut, bins = 100)
        self.le_hist = le_hist
        self.electron_cut = electron_cut
        
        return (self.le_hist, self.electron_cut)
    






