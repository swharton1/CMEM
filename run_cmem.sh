#!/bin/bash

#To compile and run this script, type this in the terminal. 
#chmod +x test.sh (only once) 
#./test.sh (this runs this script) 

echo Shell script is working...
#sleep 1

#Now try to run the python script from here. This will print all the print statements in the python code. 
#python3 test.py Martin 

#sleep 1

#To direct all the standard output (print statements) to a variable, do. However, all print statements will be combined as one string.  
#result=$(python3 test.py Geraldine) 

#echo $result 

#Let's see if I can run the run_cmem.py file. 
#echo Trying run_cmem.py 

#First, you will need to set the following environment variables for CMEM to work. You could put these in your .bashrc file and then type 
#source ~/.bashrc
export PLOT_PATH="$HOME/Code/plots/CMEM_plots/"
export PICKLE_PATH="$HOME/Code/pickled_files/CMEM_pickled_models/"
export PPMLR_PATH="/data/sol-ionosphere/SMILE/PPMLR/"

#These are the arguments for the read_ppmlr file. 
filename="S05D05V400B0000-05rad.fits"
xmin=-5
xmax=25
ymin=-25
ymax=25
zmin=-25
zmax=25


#This is the name of the emissivity model you want. jorg or cmem. 
model=cmem

#This is the name of the cost function you want to use when fitting. 
cost_func=normalised

#This is the method you want to use from the CMEM paper to initialise the starting parameters. 
init_method=1

#This is the name of a pickled model file you would read if you wanted to plot the fitted output. 
pickled_file=$filename'_'$model'_'$cost_func'_'$init_method'.pkl'

echo $pickled_file

python3 run_cmem.py $filename $xmin $xmax $ymin $ymax $zmin $zmax $model $cost_func $init_method $pickled_file

#Use the linux 'display' command to display any plots you produce. 


