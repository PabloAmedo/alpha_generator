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
from bragg_peak import*
from Clusters_calculation import*


#Lets define the classes
class Source:
    
    """This is the source class. It includes all the relevant information about the source"""
    
    def __init__(self,rate=500,energy=5.5,radius=0.35,M=933,range_alpha=4.8, We=25.7e-6):
        
        self.rate=rate #In Hz
        
        self.energy=5.5 #In Mv
        
        self.radius=radius # In cm
        
        self.range_alpha=range_alpha
        
        self.M=M
        self.z=2
        self.n_e=int(energy/We)                                  #Add We at some point
    
    def produce_alpha(self,n,store,ionization_profile,phi_in=None,ath_in=None,theta=None, n_electrons=214007): 
        """This method produces an alpha track from the alpha_tracks class"""        
        
        if ionization_profile=='Bragg':
            self.x,self.Sp,self.acum=bragg_peak(self.energy,self.range_alpha)
        else:
            self.x=None
            self.Sp=None
            self.acum=None

        for i in range(n):
            if phi_in==None:
                phi=np.random.rand()*2*np.pi
            else:
                phi=phi_in
                
            if ath_in==None:
                ath=np.random.rand()*np.pi/2
            else:
                ath=ath_in
                            
            theta=np.random.rand()*2*np.pi

            x0=np.cos(theta)*self.radius*np.random.rand()
            y0=np.sin(theta)*self.radius*np.random.rand()
                 
            alpha=Alphas_tracks(x=self.x,acum=self.acum,range_alpha=self.range_alpha,phi=phi,ath=ath,x0=x0,y0=y0,ionization_profile=ionization_profile, n_electrons=n_electrons)
            store.append(alpha)
            
        return store
    
    def collimator(self, stored_tracks, collimated_tracks, phi_col, slit, dist):
        
        """
        Method to set a collimation for the alphas. 
        
        stored_tracks       --> list with tracks previously generated
        collimated_tracks   --> list where the collimated tracks will be stored
        phi                 --> angle where the collimator is set
        slit                --> distance of the window
        dist                --> distance from the alpha to the window
        """
        
        # beta. Angle between source and window's edge
        
        beta=np.arctan((slit/2/dist))
        
        for track in stored_tracks:
            
            if track.phi<=phi_col+beta and track.phi >= phi_col-beta:
                
                collimated_tracks.append(track)
        
        return collimated_tracks



class muon_generator:
    
    """
    muon_generator class will generate 'n' muon tracks 
    with the corresponding number of ionization clusters + e-. Tracks are generated 
    randomly in a choosen rectanglar volume. We make the next assumptions:
        -We are considering that muon flux is coming from one lateral only.
        -We are only considering muons that go throught the entire region.
        -Not considering muons generated inside due to interactions.
        -Not considering muons stopping in the chamber region.
        -Muons do not loss energy
    """
    
    def __init__(self, energy=5e3,  geometry = [1,1,1]):    # Valores preliminares 
        """
        We set the general values for every muon we would like to generate.
        (This can be edited by calling the class with speciffic info about the params)
        
        example:
            muon_generator(rate=2, energy=5e3, dEdx=2.0)
        
        Sources for the used data:
        
        energy      -->         4 GeV according to (https://pdg.lbl.gov/2019/reviews/rpp2019-rev-cosmic-rays.pdf) preliminar value
        *dEdx        -->         2.061 MeV·cm2/g according to (https://pdg.lbl.gov/2017/AtomicNuclearProperties/MUE/muE_argon_gas_Ar.pdf)
                                (parametrizar como una recta dp del material y energia ???????)
        geometry    -->         chamber dimensions (maybe this fits better in other class)
        *n_cl_cm     -->         based on muon_ionization results
        
        
        *We are calculating them directly in the next step as they depend on the energy.
        """
        
        self.energy = energy  #GeV
        self.xmax, self.ymax, self.zmax = geometry
        
        
        
    def produce_muon(self,n, store, y0_in=None, phi0_in=None, theta0_in=None, e_cut = 205, gas = 'Argon', line = False):#FIXME: I've set the 104 max e cut off based on Arogancia09 paper. My calculations said 125 (kinda the same)
        
        """
        This method is used for generate n muon tracks by randomly generate a 
        y0 position where the muon will enter the chamber and also phi0 and 
        theta0 defining the direction the muon will follow across the chamber.
        
        Some considerations:
            
            -We are considering that muon flux is coming from one lateral only.
            -We are only considering muons that go throught the entire region.
            -Not considering muons generated inside due to interactions.
            -Not considering muons stopping in the chamber region.
            
            FIXME
        ^^^                                                   ^^^
        ||| (Those could be things to implement in a future)  |||
        
        -----------------------------------------------------------------------
        
        You can set:
            
            -n                 : (int) number of tracks to generate
            -store             : (list / array-like) to store tracks
            -y0                : (float) 
            -phi0              : (float)  initial position
            -theta0            : (float) 
            -e_cut             : (int) max number of e- produced in a single cl
        -----------------------------------------------------------------------
        
        Output: This method will call another one 'muon_tracks' where all the 
        relevant info will be stored (initial position, final position, n_elec,
                                      n_cl, more?)
        
        """
        if gas == 'Argon':
            gammas = np.array([4.0,3.5])
            cl_meas = np.array([27.8,28.6])
            Wi = 26.4e-6
            
        elif gas == 'Xenon':
            gammas = np.array((4.0))
            cl_meas = np.array((44))
            Wi = 21.7e-6
        dEdx, n_cl_cm = clusters_cm(Emuon = self.energy, Wi = Wi, gas = gas, gammas=gammas, cl_measurements = cl_meas )
        n_cl_cm = 34.8 #FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        track_len=[]  
        
        for i in range(n):
            #Generate the initial positions
            #FIXME
            #CHANGING THINGS TO CALCULATE THE RESOLUTION (SAME LEN)
            y0 = np.random.rand() * self.ymax
            theta0 = np.random.rand() * np.pi
# =============================================================================
            if line == True:
                y0 = self.ymax / 2
                #theta = np.pi/2
# =============================================================================
            #phi0 = np.random.rand() * 2 * np.pi 
            phi0 = 0                    #just for testing (simple case)
            """
            FIXME
            
            We are going to assume that the tracks are lines, so we can compute 
            their lengths and also their ionization using the paramentes 
            previosly obtained.
            
            
                    X
            -------------------
            |                 |
            |                 | Y
            |                 | 
            |                 |
            -------------------
            
            """
            
            
            if theta0 > np.pi/2:
                alpha = np.pi - np.pi/2 - (np.pi-theta0)
                x =  y0 / np.tan(alpha)
                
            if theta0 <= np.pi/2:
                alpha = np.pi - np.pi/2 - theta0
                x = (self.ymax - y0) / np.tan(alpha)
                
            #I should add condition about theta0 here too
            if x <= self.xmax and theta0 <= np.pi/2:
                yout = self.ymax
                xout = x
            elif x <= self.xmax and theta0 > np.pi/2:
                yout = 0
                xout = x
            elif x > self.xmax and theta0 <= np.pi/2:
                yout = self.xmax * np.tan(alpha)
                xout = self.xmax
            elif x > self.xmax and theta0 > np.pi/2:
                alpha = theta0 - np.pi / 2
                yout = self.xmax * np.tan(alpha)
                xout = self.xmax
            
            if line == True:
                yout = y0
                xout = self.xmax
            tr_len = np.sqrt(xout**2 + (yout-y0)**2)
            #Calculate number of electrons generated
            #Calculate the avg n of clusters for the track
            n_cl_avg = np.random.poisson(lam = tr_len * n_cl_cm)          #Poisson distr centered in the n_cl_len obtained from muon_ionization
            
                
            #Define the 1/n^2 distribution
            def distr_1n2(n):
                """
                exp_prob = np.array([65.6, 15.0, 6.4, 3.5, 2.25, 1.55, 1.05, 0.81, 0.61, 0.49, 0.39, 0.3,
                                 0.25, 0.2, 0.16, 0.12, 0.095, 0.075, 0.063]) / 100
                if n < 19:
                    return exp_prob[n]
                if n >= 19:
                    return 21.6/n**2"""
                return 1/n**2

            n_e_cl=[]                   #number of e- generated in each cluster
            """
            FIXME
            
            The next way to compute the distribution could be very slow
            """
            while len(n_e_cl) < n_cl_avg:
                a = np.random.randint(1, e_cut)
                b = np.random.rand()
                
                if b < distr_1n2(a):
                    n_e_cl.append(a)
                    
            #Calculation of initial and final position
            """ 
            FIXME
            This is just true for the 2D case, change for more generality
            """
            
            z0 = tr_len * np.sin(theta0) * np.cos(phi0)
            y0 = y0
            x0 = 0
            
            z = z0
            y = yout
            x =xout
            
            #Uniform distribution of the clusters along the track
            cl_position = np.random.uniform(low = 0, high = xout, size = len(n_e_cl))
            cl_position = np.sort(cl_position)
            
            
            #saving length
            muon = muon_tracks(energy = self.energy, track_length = tr_len, initial_position = [x0, y0, z0], 
                               final_position =[x, y, z] , clusters = cl_position, n_e_cl = n_e_cl, geometry = [self.xmax, self.ymax, self.zmax])
            store.append(muon)
            
        return n_e_cl, cl_position


class muon_tracks(muon_generator):
    
    """
    muon_tracks class. It contains all the relevant information about the muons 
    travel across the chamber (initial position, final position, length, 
                               n_cl, n_e, ...)
    """
    
    def __init__(self, energy, track_length, initial_position, final_position, clusters, n_e_cl, geometry = [1,1,1], spread = 0):
        
        self.energy = energy
        self.track_length = track_length
        self.x0, self.y0, self.z0 = initial_position
        self.x, self.y, self.z = final_position
        self.clusters = clusters    #z coord of each cluster
        self.n_e_cl = n_e_cl
        self.n_electrons = np.sum(n_e_cl)
        self.xmax, self.ymax, self.zmax = geometry
        
        if self.y > self.y0:
            self.phi = np.arctan(abs(self.y-self.y0) / (self.x-self.x0))
        elif self.y == self.y0:
            self.phi = np.pi / 2
        else:
            self.phi = np.pi / 2 + np.arctan(abs(self.y-self.y0) / (self.x-self.x0))
        
        #This is defined at 2 sigma and will be handled by the Difussion_handler class
        self.spread=spread
        
        #Storage of electrons (z,y) positions
        self.electron_positions=np.zeros([self.n_electrons,2]) # coordenadas [x,y] para los 50 electrones
        
        #Storage of electrons (z,y) positions after diffusion
        self.electron_positions_diff=np.zeros([self.n_electrons,2])
        
    def fill(self):
        """
        This method is used to fill the clean muon track through the chamber 
        with the electrons produced in the ionization.
        """
        
        #Calculate the slope (assuming 2D trakcs)
        #y(x) = m * x + n
        slope = (self.y-self.y0) / (self.x-self.x0)
        
        def lineal_track(x):
            return slope * x + self.y0
        
        cl_y_tr = []                  #List to store y coord of cl around track 
        cl_x_tr = []                  #List to store x coord of cl around track

        for i in range(len(self.clusters)):
            cl_y_pos = lineal_track(self.clusters[i])
            if cl_y_pos <= self.ymax and cl_y_pos >= 0:
                cl_y_tr.append(cl_y_pos)
                cl_x_tr.append(self.clusters[i])
        
        #Assign positions to e in cluster
        aux_electrons_total = 0
        for i in range(len(self.clusters)):
            for j in range(self.n_e_cl[i]):
                #Change the positions
                self.electron_positions[aux_electrons_total,0] = cl_x_tr[i]
                self.electron_positions[aux_electrons_total,1] = cl_y_tr[i]
                #Change the positions
                self.electron_positions_diff[aux_electrons_total,0] = cl_x_tr[i]
                self.electron_positions_diff[aux_electrons_total,1] = cl_y_tr[i]
                
        
                # Calculate the transverse diffusion on the z-y plane
                if self.spread!= 0:
                    #Now we need to draw the number from the gaussian distribution
                    pos = np.random.multivariate_normal((0,0), cov=self.spread*np.identity(2),size=1)
                    
                    self.electron_positions_diff[aux_electrons_total,0]+=pos[:,0]
                    
                    self.electron_positions_diff[aux_electrons_total,1]+=pos[:,1]
                    
                aux_electrons_total = aux_electrons_total + 1
        
        return cl_x_tr, cl_y_tr
    
    def select(self,variable_scan="x",min_var=0,max_var=100,array_to_store=None):
        """This method scans the positions of the electrons either in the x or y direction and
        it selects them based on a condition on the other variable. Ex: we scan in x and select
        electrons whose y position fulfills y_min<y<y_max"""
        
        #Select the apropiate index, the opossite of variable_scan
        if variable_scan=="x":
            index=1;other_index=0
        else:
            index=0;other_index=1
            
        #Iterate over the electrons
        for j in range(self.n_electrons):
            #Check if they are whitin our current 
            if min_var<self.electron_positions_diff[j,index]<=max_var:
                array_to_store.append(self.electron_positions_diff[j,other_index])
                
        return array_to_store
    
        
class Alphas_tracks(Source):
    
    """This is the alpha class. It contains all the relevant information about the alpha track,
    like the ionization profile, positions of electrons, etc"""
    
    def __init__(self,x,acum,range_alpha,phi=None,ath=None,x0=None,y0=None,spread=0,ionization_profile="Flat", n_electrons=int(5.5/25.7e-6)):
        
        self.ionization=ionization_profile
        self.X=x
        self.acum=acum            #si no da error
                                  #son vectores x y Sp, necesarios para los randoms
        
        #This is the range of the alpha. Get it from NIST
        self.range_max=range_alpha #In cm
        
        #Athimutal angle. This reduces our effective alpha range by the projection
        self.ath=ath
        
        #Define an effective range
        self.range=self.range_max*np.cos(ath)
                      
        self.x0=x0
        self.y0=y0
      
        #X position in the plane
        self.x= np.cos(phi)*self.range*np.cos(ath) + self.x0 
        
        #Y position in the plane
        self.y=np.sin(phi)*self.range*np.cos(ath)  + self.y0
        
        #Phi angle
        self.phi=phi
               
        #This is defined at 2 sigma and will be handled by the Difussion_handler class
        self.spread=spread
        
        #Set the ionization profile to be flat between the 0 and the range
        self.dict_ion_prof={
            "Flat":np.random.rand,
            "Bragg":random_bragg
            }
        
        self.ionization_profile=self.dict_ion_prof[self.ionization]
        
        #Number of electrons in a track
        # self.n_electrons=source.energy/Gas.I
        # self.n_electrons=50
        self.n_electrons=n_electrons
        
        #Storage of electrons (x,y) positions
        self.electron_positions=np.zeros([n_electrons,2]) # coordenadas [x,y] para los 50 electrones
        
        #Storage of electrons (x,y) positions after diffusion
        self.electron_positions_diff=np.zeros([n_electrons,2])
        
    def fill(self):
        
        """This method must only be called when we want to fill the final alpha track with
        electrons. 
        
        First it picks a random number following a given ionization profile to situate the electron
        along the track. Then it picks a number to situate along transverse gaussian distribution
        that mimics diffusion. If the spread is zero then we don't apply the transversal diffusion.
        """
        
        #Draw a number from the ionization profile distribution
        #self.ionization_profile=1
        for i in range(self.n_electrons):
            
            #Get the radial position of the electron alongside the track
            if self.ionization=='Bragg':
                r_pos=self.ionization_profile(self.X,self.acum)*self.range
            else:
                r_pos=self.ionization_profile()*self.range
                            
            #Change the positions
            self.electron_positions[i,0]=np.cos(self.phi)*r_pos + self.x0
            self.electron_positions[i,1]=np.sin(self.phi)*r_pos + self.y0
            #Change the positions
            self.electron_positions_diff[i,0]=np.cos(self.phi)*r_pos + self.x0
            self.electron_positions_diff[i,1]=np.sin(self.phi)*r_pos + self.y0
            
            # Calculate the transverse diffusion on the x-y plane
            if self.spread!=0:  
                #Now we need to draw the number from the gaussian distribution
                pos = np.random.multivariate_normal((0,0), cov=self.spread*np.identity(2),size=1)
                
                self.electron_positions_diff[i,0]+=pos[:,0]
                
                self.electron_positions_diff[i,1]+=pos[:,1]
                
                
    def select(self,variable_scan="x",min_var=0,max_var=100,array_to_store=None):
        """This method scans the positions of the electrons either in the x or y direction and
        it selects them based on a condition on the other variable. Ex: we scan in x and select
        electrons whose y position fulfills y_min<y<y_max"""
        
        #Select the apropiate index, the opossite of variable_scan
        if variable_scan=="x":
            index=1;other_index=0
        else:
            index=0;other_index=1
            
        #Iterate over the electrons
        for j in range(self.n_electrons):
            #Check if they are whitin our current 
            if min_var<self.electron_positions_diff[j,index]<=max_var:
                array_to_store.append(self.electron_positions_diff[j,other_index])
                
        return array_to_store
    

     
class Gas:
    
    """This is the gass class. It includes information about the W, drift velocity, etc"""
    
    def __init__(self,vz=None,Dt=None,Dl=None,Wi=None,density=None):
        
        self.vz=vz
        self.Dt=Dt
        self.Dl=Dl
        self.Wi=Wi
        self.density=density
        
        
        self.A=40
        self.Z=18
        self.I=25.7e-6 #Mev

class Diffusion_handler(Gas):

    """This class handles the application of diffusion from the gas. It inherits from the gass
    class"""
    def __init__(self,sigma_diff=0.25,sigma_PSF=0):
        
        self.u = 1
        self.sigma_diff = sigma_diff
        self.sigma_PSF = sigma_PSF
        
    def diffuse(self,alpha_list):
        
        #Take an alpha track and add a difussion in cm
        for track in alpha_list:
            track.spread= ( self.sigma_diff**2 + self.sigma_PSF )**0.5  # Desviación estándar de la Gaussiana (7.5 mm??? too much???)
    


    
class Noise():

    """This class handles the noise addition to the final 2D histogram where we assume we have
    our array of pixels
    
    We care about the spread of the noise, not so much about the absolute value, so this noise
    magnitudes will probably be about spread
    
    """
    
    def __init__(self,dark_current=190):
        
        #Dark_current (electrons/s) of the camera is 190 e/s. This depends exponentially on the
        #temperature
        self.dark_current=dark_current
    
    def add_noise(self,exposition_time,Hist2d,sigma=5.7):
        
        """This method adds noise to each pixel of a 2D histogram as a function of the
        electronic noise and the exposure time"""
        
        #Electronic noise (e) = dark_current (e/s) * exposition_time (s)
        #These units are in electrons
        electronic_noise=self.dark_current*exposition_time
        
        #Create an array of n_bins*n_bins and populate them with random numbers
        #random_gauss=np.random.normal(scale=electronic_noise,size=Hist2d.shape)
        random_gauss=np.random.normal(loc=electronic_noise,size=Hist2d.shape,scale=sigma)
        #Ignore the negative ones?
        random_gauss=abs(random_gauss)
        
        #Update the 2d histogram with this values and return it
        
        return Hist2d+random_gauss
    
    def add_noise_abs(self,Hist2d,mean=0,sigma=5.7):
        
        """This method adds noise to each pixel of a 2D histogram using a gaussian distribution
        for each pixel"""
        
        random_gauss=np.random.normal(loc=mean,size=Hist2d.shape,scale=sigma)
        #Ignore the negative ones?
        random_gauss=abs(random_gauss)
        
        #Update the 2d histogram with this values and return it
        
        return Hist2d+random_gauss
        
        

class Image_2D():
    
    """This class is essentially a 2D histogram that also stores extra information at the truth level:
        how many tracks were produced, their positions, the positions of the electrons, etc. This
        is done with numpy.histogram2d and the axis are already inverted..
        
        This object will be exported and loaded by our analysis framework"""
        
    def __init__(self,track_list,hist_args={"bins":10},QE=1,GE=1,Tp=1,gain=1,ph_poiss=False):
        
        #List of tracks
        self.track_list=track_list
        
        #Get a list of all electron positions in all tracks
        #self.x_pos = np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,0] for i in range(len(track_list))]))
        #self.y_pos = np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,1] for i in range(len(track_list))]))
        
        self.x_pos = np.concatenate([track_list[i].electron_positions_diff[:,0] for i in range(len(track_list))])
        self.y_pos = np.concatenate([track_list[i].electron_positions_diff[:,1] for i in range(len(track_list))])
        
        
        #Get a list of all the true electron's positions in all tracks
        #self.x_pos_true=np.ndarray.flatten(np.asarray([track_list[i].electron_positions[:,0] for i in range(len(track_list))]))
        #self.y_pos_true=np.ndarray.flatten(np.asarray([track_list[i].electron_positions[:,1] for i in range(len(track_list))]))
        
        self.x_pos_true = np.concatenate([track_list[i].electron_positions[:,0] for i in range(len(track_list))])
        self.y_pos_true = np.concatenate([track_list[i].electron_positions[:,1] for i in range(len(track_list))])
        
        
        
        #Define the quantum efficiency, the geometrical efficiency, the transparency and the gain
        self.QE=QE
        self.GE=GE
        self.Tp=Tp
        self.gain=gain    
    
        #Get the 2d hist for all the electrons and the edges from the list of tracks.
        #Extra arguments can be passed to the 2D numpy histogram function
        self.Hist2D_e,self.x_edges,self.y_edges=np.histogram2d(x=self.x_pos,y=self.y_pos,**hist_args)
        #Since numpy inverts the (y,x) histogram, invert it again
        self.Hist2D_e=self.Hist2D_e.T
        
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
            
    def plot_hist(self,fig_in=None,axis_list_in=None,noise_object=None,exposition_time=0):
        """This simply plots the 2D histogram"""
    
        H,yedges,xedges=self.Hist2D,self.y_edges,self.x_edges
        
        #If a figure and axis are provided, use them. Otherwise come up with our own
        if fig_in==None or axis_list_in==None:
            #Create a figure
            fig, (ax1,ax2) = plt.subplots(ncols=2)
            
        else:
            fig=fig_in
            ax1=axis_list_in[0]
            ax2=axis_list_in[1]
            
            
        #Plot it
        ax1.pcolormesh(xedges, yedges, H, cmap='rainbow',vmin = np.min(H), vmax = np.max(H))
        
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
    
    def plot_x(self):
        electron_cut=[];variable_cut='x';other_variable='y'
        for track in self.track_list:
            electron_cut=track.select(variable_cut,min_var=-0.5,max_var=0.5,array_to_store=electron_cut)
        
        le_hist=np.histogram(electron_cut,bins=100)
        self.le_hist=le_hist
        self.electron_cut=electron_cut
        
        return (self.le_hist,self.electron_cut)
    
    
    
