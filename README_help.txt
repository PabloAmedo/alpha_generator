QUICK DOCUMENTATION TO RUN THE SIMULATIONS:


Things that you can edit and things you should not :)


***WARNING***
Some scripts could raise errors because I moved them to other folders!!



FILES:


-The 'example_...py' files are examples of how to run the main scripts with the classes and functions deffinitions. It's where you can define your parameters and values for the simulation. 		--->	 I'd recommend to copy them and work from there to keep a useful example alive. 

-Other '.py' files are the scripts where the classes and funtions are defined. It has some default values for the parameters.
It's recommended to do not touch them if you only want to run a simulation, pls!!  :) 
Just consult the parameters that you have to change and edit them in the examples files!!!

-The .csv files are files with stored data that it's being called by the scripts. Mostly they contain info about the camara that has been used to take the data (datos.csv), so if you want to change that info you can rewrite the file with the desire info.

	cal_hole,cal_radius,cal_fab ---> those are the ones which needs a review !!!!!

-The .tiff files are the images taken with the camera so they are the 'real' data to compare/analyze. You will need to specify the path + file in your script (by default it is called in the example_optical_gain.py file). 



----------------------------------------------------------------------------------------------------------------------------



A more detailed documentation will be added to this folder with speciffic info about classes and methods in a future.