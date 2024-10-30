print ("Running CMEM Code")

#Import required python modules here. 
import sys 
import os 
import CMEM

#Argument unpacking is probably going to be the first step. These arguments will either be set in shell (for now), or later in IDL. 
script = sys.argv[0]
filename = str(sys.argv[1]) 
xmin = float(sys.argv[2])
xmax = float(sys.argv[3])
ymin = float(sys.argv[4])
ymax = float(sys.argv[5])
zmin = float(sys.argv[6])
zmax = float(sys.argv[7])
model = str(sys.argv[8])
cost_func = str(sys.argv[9]) 
init_method = int(sys.argv[10])
pickle_file = str(sys.argv[11]) 

#Check that the environment variables have been set and actually exist.
PPMLR_PATH = os.environ.get("PPMLR_PATH")
PLOT_PATH = os.environ.get("PLOT_PATH") 
PICKLE_PATH = os.environ.get("PICKLE_PATH") 

if not os.path.isdir(PPMLR_PATH): 
	raise ValueError('PPMLR_PATH={} does not exist. Set the variable correctly.'.format(PPMLR_PATH))
if not os.path.isdir(PLOT_PATH): 
	raise ValueError('PLOT_PATH={} does not exist. Set the variable correctly.'.format(PLOT_PATH))
if not os.path.isdir(PICKLE_PATH): 
	raise ValueError('PICKLE_PATH={} does not exist. Set the variable correctly.'.format(PICKLE_PATH))	


#Set arguments here. These may come from IDL or Shell so can be commented out when needed. 
#filename="S05D05V400B0000-05rad.fits"
#xmin=-5
#xmax=25
#ymin=-25
#ymax=25
#zmin=-25
#zmax=25
#model='cmem'
#cost_func='normalised'
#pickle_file=''


#params is None by default and it uses method 2 to initialise parameters. You can set them manually here. 
params0=None


#READING IN THE PPMLR DATA CUBES. 
#################################

#This will read in the PPMLR emissivity cube. 
print ("Filename is: ", filename) 
print ("Get PPMLR data") 
ppmlr = CMEM.ppmlr_fits.read_ppmlr_fits(filename=filename, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
print ("Density = ", ppmlr.density) 

#PLOTTING THE PPMLR DATA CUBES.
###############################

#This will plot the XY and XZ planes. 
#ppmlr.plot_both_planes(cmap="hot", levels=100, vmin=-8, vmax=-4, save=True, savetag="")

#PLOTTING A NON-OPTIMISED MODEL AGAINST THE PPMLR CUBES.
########################################################

#This code is how to read the code for plotting the non-optimised models against the PPMLR simulation. 
#There are two ways to run it. You can pass the filename in and it will read the PPMLR file, or you can pass in the PPMLR object you might have already created above. 

#Create cdm object. 
#cdm = CMEM.visualise_nonopt.compare_data_model(filename=filename, ppmlr=ppmlr, model=model, params0=params0, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)

#How to create a plot showing 2D planes through the simulation and the model with the initial parameters selected. 
#cdm.plot_planes(cmap='hot', vmin=-8, vmax=-4, levels=100, save=True, savetag="")

#How to create a plot showing the values along the Earth-Sun line. 
#cdm.plot_earth_sun_line(save=True, savetag="") 

#FIT A MODEL TO THE PPMLR DATA CUBES. 
#####################################
print ('Fitting...') 
#This code is for fitting a model to the PPMLR simulation. 
fit_model = CMEM.fit_emissivity_models.threed_models(filename=filename, ppmlr=ppmlr, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
#print ("Start fitting...")
#This does the fitting. 
fit_model.fit_function_with_nelder_mead(model=model, params0=params0, set_param_bounds=False, cost_func=cost_func, init_method=init_method)

#Save the model output to a pickle file. 
fit_model.write_pickle(fname=pickle_file, savetag="") 

#PLOT THE OUTPUT OF FITTING A MODEL TO THE PPMLR DATA CUBES. 
############################################################

#Visualisation tools for the fitted output. This reads the pickled file. 
#analysis = CMEM.visualise_models.analyse_model(filename=pickle_file)

#To plot the parameter variation with iteration number. 
#analysis.plot_change_in_parameters(save=False, savetag="") 

#To plot the optimised model against the simulation. 
#analysis.plot_planes(cmap='hot', vmin=-8, vmax=-4, levels=100, save=True, savetag="")

#To plot along the sun-earth line. 
#analysis.plot_earth_sun_line(save=True, savetag="") 

#THERE ARE SOME FURTHER ANALYSIS FUNCTIONS IN COMPARE_OPTIMISED_MODELS BUT THESE REQUIRE ALL SIMULATIONS TO BE RUN. RECOMMEND DOING THIS YOUR OWN WAY. 

#These produce plots comparing minimum cost values between models. 
#com = CMEM.compare_optimised_models.compare_models()
#com.compare_cost_between_models(save=True, savetag="", set_A1A2=False, add_bad_fits=True)
#com.compare_cost_between_models_both_methods(save=True, savetag="", add_bad_fits=False)

#These produce plots the relationships with SW density for the first four parameters.  
#pr = CMEM.compare_optimised_models.parameter_relationships(model=model)
#pr.plot_A1_A2_relationships(save=True, savetag="")

#This produces the plot for the relationship between SW density and all optimised parameters. 
#opr = CMEM.compare_optimised_models.optimal_parameter_relationships(model='cmem')
#opr.plot_optimal_parameters(save=True) 

#This produces the plot showing all the definitions of the magnetopause boundary. 
#mag = CMEM.compare_optimised_models.magnetopause_model()
#mag.plot_earth_sun_line(sim=3, xlim=[7,10], save=True)
#mag.plot_all_boundaries(save=True)



