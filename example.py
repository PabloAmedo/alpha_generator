# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 10:48:14 2023

@author: jacob
"""

from alpha_track_simulator import*



#==This is the main program==
#Create an argon object from the gas class
argon=Gas(1,1,1,1,1)
#Create a source and a list to store

source=Source(radius=0.1)
track_list=[]

#Creat some alpha tracks and set some parameters
n_tracks=10;ath_angle=0

exposition_time=n_tracks/source.rate #In seconds


source.produce_alpha(n=n_tracks,store=track_list,ath_in=None,ionization_profile="Bragg")


#Create a difussion handler and change the diffusion of the tracks
diff_handler=Diffusion_handler()
diff_handler.diffuse(track_list)

#Generate the electrons alongside the tracks
for i in track_list: i.fill()


#Create a noise object with a given dark noise
noise=Noise(50)

#Plot the tracks
image2d=Image_2D(track_list=track_list,hist_args={"bins":20})
# #Plot the tracks
image2d.track_plot()
# image2d.plot_hist(noise_object=noise,exposition_time=exposition_time)




















