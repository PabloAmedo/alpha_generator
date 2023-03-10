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
        electrons. If the spread is zero then we don't apply the transversal diffusion
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
    
    
    
#==This is the main program==
#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)
#Create a source and a list to store

# ______Jacobo_________________________________________________________________
source=Source(radius=0.1)
track_list=[]
#Creat some alpha tracks
n_tracks=10;ath_angle=0
source.produce_alpha(n=n_tracks,ath_in=None,store=track_list)
#______________________________________________________________________________

#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler()
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()


#===Plot section===
#This is super slow if you have to do it track by track. Just skip it if there are too many
if n_tracks<100:
    #Plot the current resut
    fig=plt.figure(figsize=(6,6))
    #fig.subplots_adjust(wspace=0.5,hspace=0.3)
    ax=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)
    #Plot the tracks
    for i in range(len(track_list)):
        
        ax.plot([track_list[i].x0,track_list[i].x],[track_list[i].y0,track_list[i].y])        

    ax.set_xlabel("x (cm) ")
    ax.set_ylabel("y (cm) ")

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

#Get a list of electron positions
x_pos=np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,0] for i in range(len(track_list))]))
y_pos=np.ndarray.flatten(np.asarray([track_list[i].electron_positions_diff[:,1] for i in range(len(track_list))]))

#Get the 2d hist
H, yedges, xedges = np.histogram2d(y_pos, x_pos, bins=30)
#Create a figure
fig, ax1 = plt.subplots(ncols=1)
#Plot it
ax1.pcolormesh(xedges, yedges, H, cmap='rainbow')

#ax1.plot(x, 2*np.log(x), 'k-')

ax1.set_xlim(x_pos.min(), x_pos.max())
ax1.set_ylim(y_pos.min(), y_pos.max())
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Nº of tracks '+ str(n_tracks)+" , equivalent to "+str(n_tracks/500)+" s exposure time")

#Try to get the nº of electrons in a cut
electron_cut=[];variable_cut="x";other_variable="y"
for track in track_list:
    #Get the cut
    electron_cut=track.select(variable_cut,min_var=-0.5,max_var=0.5,array_to_store=electron_cut)
    
#Create a figure
le_hist=np.histogram(electron_cut,bins=15)
fig, ax1 = plt.subplots(ncols=1)
ax1.hist(electron_cut,bins=15,fill=False,label="Montecarlo")
ax1.set_xlabel("x (cm)")
ax1.set_ylabel("N")

#Compare with some of Jacobo's data
path="C:/Users/jacob/OneDrive - Universidade de Santiago de Compostela/Documentos (Escritorio)/Física/Laboratorio/Imágenes Teledyne/Fat gem 8 to 0 medidas/corte x/" 
data_df=pd.read_csv(path+"sec x data.csv")
# Hay 58 px entre agujeros y 5 mm entre agujeros
cal = 58/5 *10 #11.6 px/mm

#ax1.bar(data_df["px"]*2/2500,data_df["1"]*max(le_hist[0])/max(data_df["1"]),fill=False)
ax1.scatter(data_df["px"]*2/2500-1.12,data_df["19"]*max(le_hist[0])/max(data_df["19"]),c="k",label="Data")
ax1.legend()
ax1.set_title(variable_cut+"-axis cut with "+str(-0.5)+"<y< "+str(0.5))
