;README
;This IDL script will execute run_cmem.py. By default all plotting functions are turned off. 
;The script does not return anything. Figures or model files will be outputted at the 
;environment variables you set below. 


;First, you will need to set the following environment variables for CMEM to work. You could put these in your .bashrc file and then type source ~/.bashrc
HOME=getenv('HOME')

;This is where the ppmlr emissivity files are located. 
setenv, 'PPMLR_PATH=/data/sol-ionosphere/SMILE/PPMLR/'

;This is where you want any plots to appear. You may need to make this directory the first time. 
setenv, 'PLOT_PATH='+HOME+'/Code/plots/CMEM_plots/'

;This is where the output of the fitting process goes, as a pickle file. 
setenv, 'PICKLE_PATH='+HOME+'/Code/pickled_files/CMEM_pickled_models/' 

;These are the arguments for the read_ppmlr file. 
filename="S05D12.25V400B0000-05rad.dat"
xmin=-5
xmax=25
ymin=-25
ymax=25
zmin=-25
zmax=25

;This is the name of the emissivity model you want. jorg or cmem. 
model="cmem"

;This is the name of the cost function you want to use when fitting. 
cost_func="normalised"

;This is the method you want to use from the CMEM paper to initialise the starting parameters. 
init_method=2

;This is the name of a pickled model file you would read if you wanted to plot the fitted output. You do not need to set this, but leave this default name in so you don't break it!  
pickled_file="S05D7.5V400B0000-05rad.dat_cmem_normalised.pkl"

;Convert any integers or floats to trimmed strings. 
xmin = strtrim(string(xmin), 1) 
xmax = strtrim(string(xmax), 1) 
ymin = strtrim(string(ymin), 1) 
ymax = strtrim(string(ymax), 1) 
zmin = strtrim(string(zmin), 1) 
zmax = strtrim(string(zmax), 1) 
init_method = strtrim(string(init_method), 1) 


;Form the command to go into the spawn command. 
cmd = strjoin(['python3 run_cmem.py', filename, xmin, xmax, ymin, ymax, zmin, zmax, model, cost_func, init_method, pickled_file], ' ')
print, cmd
spawn, cmd 
