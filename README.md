# CMEM
This is the code for the Cusp Magnetosheath Emissivity Model (CMEM) 

All figures and analysis in the CMEM paper can be created with this code. 

The code can be run from the shell, ipython or IDL. 

All the code is in the CMEM folder. 

To learn how to run the code, I recommend you look in run_cmem.py. It contains commented out examples of how to run everything in there. You can either type the commands into ipython3 or create your own version of run_cmem.py to just run the commands you want. 

THINGS YOU NEED TO DO: 

1. Install SXI_Core first from my github site. Follow those instructions. This contains some core, common functions across all of my projects. 

2. You may want to add CMEM to your PYTHONPATH variable, say in your .bashrc file. Then you can call it from anywhere. Mine looks like:
PYTHONPATH=$PYTHONPATH:~/Code/CMEM/

3. There are a series of paths you will need to set in CMEM/__init__.py. read_fits_cube uses paths stored as environment variables to find the emissivity cubes. You need to set these in the init file. Change whichever of these variables you need too: PPMLR_PATH, OPENGGCM_PATH or BATSRUS_PATH. You will also need to set the default location for your PLOT_PATH and PICKLE_PATH variables, where the output goes. 
