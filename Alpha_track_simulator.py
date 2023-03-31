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


#Lets define the classes
class Source:
    
    """This is the source class. It includes all the relevant information about the source"""
    
    def __init__(self,rate=500,energy=5.5,radius=1): #Jacobo
        
        self.rate=rate #In Hz
        
        self.energy=5.5 #In Mv
        
        self.radius=radius # In cm
        
    def produce_alpha(self,n,store,phi_in=None,ath_in=None,theta=None): #Jacobo
        """This method produces an alpha track from the alpha_tracks class"""
        
        for i in range(n):
            if phi_in==None:
                phi=np.random.rand()*2*np.pi
            else:
                phi=phi_in
                
            if ath_in==None:
                ath=np.random.rand()*np.pi/2
            else:
                ath=ath_in
                
    # ______Jacobo_____________________________________________________________
            
            theta=np.random.rand()*2*np.pi

            x0=np.cos(theta)*self.radius*np.random.rand()
            y0=np.sin(theta)*self.radius*np.random.rand()
                 
            alpha=Alphas_tracks(phi=phi,ath=ath,x0=x0,y0=y0)
            store.append(alpha)
    # _________________________________________________________________________            
            
        return store
        
class Gas:
    """This is the gass class. It includes information about the W, drift velocity, etc"""
    
    def __init__(self,vz=None,Dt=None,Dl=None,Wi=None,density=None):
        
        self.vz=vz
        self.Dt=Dt
        self.Dl=Dl
        self.Wi=Wi
        self.density=density
            
class Alphas_tracks:
    
    """This is the alpha class. It contains all the relevant information about the alpha track,
    like the ionization profile, positions of electrons, etc
    
    
    """
    
    def __init__(self,range_alpha=1,phi=None,ath=None,x0=None,y0=None,spread=0,ionization_profile="Flat"): #Jacobo
                
        #This is the range of the alpha. Get it from NIST
        self.range_max=range_alpha #In cm
        
        #Athimutal angle. This reduces our effective alpha range by the projection
        self.ath=ath
        
        #Define an effective range
        self.range=self.range_max*np.cos(ath)
                
        # ______Jacobo_________________________________________________________
      
        self.x0=x0
        self.y0=y0
      
        #X position in the plane
        self.x= np.cos(phi)*self.range*np.cos(ath) + self.x0 
        # + fuente.radius*np.random.rand() 
        
        #Y position in the plane
        self.y=np.sin(phi)*self.range*np.cos(ath)  + self.y0
        # + fuente.radius*np.random.rand()  
        # _____________________________________________________________________
        
        #Phi angle
        self.phi=phi
               
        #This is defined at 2 sigma and will be handled by the Difussion_handler class
        self.spread=spread
        
        #Set the ionization profile to be flat between the 0 and the range
        self.dict_ion_prof={"Flat":np.random.rand}
        
        self.ionization_profile=self.dict_ion_prof[ionization_profile]
        
        #Number of electrons in a track
        self.n_electrons=50
        
        #Storage of electrons (x,y) positions
        self.electron_positions=np.zeros([self.n_electrons,2]) # coordenadas [x,y] para los 50 electrones
        
        #Storage of electrons (x,y) positions after diffusion
        self.electron_positions_diff=np.zeros([self.n_electrons,2])
        
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
            r_pos=self.ionization_profile()*self.range
            #Change the positions
            
            # ______Jacobo_____________________________________________________
            
            # En principio debería también de cambiar las posiciones de estos 
            # electrones
            
            self.electron_positions[i,0]=np.cos(self.phi)*r_pos + self.x0
            self.electron_positions[i,1]=np.sin(self.phi)*r_pos + self.y0
            #Change the positions
            self.electron_positions_diff[i,0]=np.cos(self.phi)*r_pos + self.x0
            self.electron_positions_diff[i,1]=np.sin(self.phi)*r_pos + self.y0
            # _________________________________________________________________

            if self.spread!=0:
                #Now we need to draw the number from the gaussian distribution
                r_pos_new=np.random.normal(loc=0,scale=self.spread)
                
                if r_pos_new>=0:
                    #Change the positions
                    self.electron_positions_diff[i,0]+=np.sin(self.phi)*abs(r_pos_new)
                    
                    self.electron_positions_diff[i,1]+=-np.cos(self.phi)*abs(r_pos_new)
                    
                elif r_pos_new<0:
                    #Change the positions
                    self.electron_positions_diff[i,0]+=-np.sin(self.phi)*abs(r_pos_new)
                    
                    self.electron_positions_diff[i,1]+=np.cos(self.phi)*abs(r_pos_new)
            

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
                

class Diffusion_handler(Gas):

    """This class handles the application of diffusion from the gas. It inherits from the gass
    class"""
    def __init__(self):

        self.u=1
        
    def diffuse(self,alpha_list):
        
        #Take an alpha track and add a difussion
        for track in alpha_list:
            track.spread=0.01
    
class Noise():

    """This class handles the noise addition to the final 2D histogram where we assume we have
    our array of pixels
    
    We care about the spread of the noise, not so much about the absolute value, so this noise
    magnitudes will probably be about spread
    
    """
    
    def __init__(self,dark_current):
        
        #Dark_current (electrons/s)
        self.dark_current=dark_current
        

    def add_noise(self,exposition_time,Hist2d):

        """This method adds noise to each pixel of a 2D histogram as a function of the
        electronic noise and the exposure time"""        

        #Electronic noise (e) = dark_current (e/s) * exposition_time (s)
        #These units are in electrons
        electronic_noise=self.dark_current*exposition_time
        
        #Create an array of n_bins*n_bins and populate them with random numbers
        random_gauss=np.random.normal(scale=electronic_noise,size=Hist2d.shape)
        
        #Ignore the negative ones?
        random_gauss=abs(random_gauss)
        
        #Update the 2d histogram with this values and return it
        
        return Hist2d+random_gauss

class Image_2D():
    
    """This class is essentially a 2D histogram that also stores extra information at the truth level:
        how many tracks were produced, their positions, the positions of the electrons, etc. This
        is done with numpy.histogram2d and the axis are already inverted..
        
        This object will be exported and loaded by our analysis framework"""
        
    def __init__(self,track_list,hist_args={"bins":10}):
        
        #List of tracks
        self.track_list=track_list

        #Get a list of all electron positions in all tracks
        self.x_pos=np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,0] for i in range(len(track_list))]))
        self.y_pos=np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,1] for i in range(len(track_list))]))
        
        #Get a list of all the true electron's positions in all tracks
        self.x_pos_true=np.ndarray.flatten(np.asarray([track_list[i].electron_positions[:,0] for i in range(len(track_list))]))
        self.y_pos_true=np.ndarray.flatten(np.asarray([track_list[i].electron_positions[:,1] for i in range(len(track_list))]))
    
    
        #Get the 2d hist and the edges from the list of tracks. Extra arguments can be passed to the
        #2D numpy histogram function
        self.Hist2D,self.x_edges,self.y_edges=np.histogram2d(x=self.x_pos,y=self.y_pos,**hist_args)
        #Since numpy inverts the (y,x) histogram, invert it again
        self.Hist2D=self.Hist2D.T
        
    
    def track_plot(self,fig_in=None,axis_list_in=None):
        
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
        if n_tracks<100:

            #Plot the tracks
            for i in range(len(self.track_list)):
                
                ax.plot([track_list[i].x0,track_list[i].x],[track_list[i].y0,track_list[i].y])        
        
            ax.set_xlabel("x (cm) ")
            ax.set_ylabel("y (cm) ")
            ax.set_title("Original tracks")
        
            #Plot the tracks
            for i in range(len(track_list)):
                
                ax2.scatter(track_list[i].electron_positions[:,0],track_list[i].electron_positions[:,1],marker="o")
        
            #Restart the color cyle of the axis
            ax2.set_prop_cycle(None)
            
            #Plot the tracks with the diffused electrons
            for i in range(len(track_list)):
                
                ax2.scatter(track_list[i].electron_positions_diff[:,0],track_list[i].electron_positions_diff[:,1],marker="x")
            
            ax2.set_xlabel("x (cm) ")
            ax2.set_ylabel("y (cm) ")
            ax2.set_title("Electron positions")
            
        return fig
            
    def plot_hist(self,fig_in=None,axis_list_in=None,noise_object=None,exposition_time=0):
        """This simply plots the 2D histogram"""
    
        H, yedges, xedges=self.Hist2D,self.y_edges,self.x_edges
        
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

#==This is the main program==
#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)
#Create a source and a list to store

# ______Jacobo_________________________________________________________________
source=Source(radius=0.1)
track_list=[]

#Creat some alpha tracks and set some parameters
n_tracks=10;ath_angle=0

exposition_time=n_tracks/source.rate #In seconds

source.produce_alpha(n=n_tracks,ath_in=None,store=track_list)
#______________________________________________________________________________

#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler()
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()

#Create a noise object with a given dark noise
noise=Noise(50)

#Plot the tracks
image2d=Image_2D(track_list=track_list,hist_args={"bins":20})
#Plot the tracks
image2d.track_plot()
image2d.plot_hist(noise_object=noise,exposition_time=exposition_time)



