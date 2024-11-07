import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import string 
from time import process_time
from matplotlib.patches import Wedge, Polygon, Circle
import pickle
import os

try: 
    #from . import read_ppmlr
    from . import ppmlr_fits
    from . import boundary_emissivity_functions as bef
    from . import set_initial_params as sip 
    from . import coord_conv as cconv 
    
except(ImportError):
    print ("Are you working from the right directory? ")
    print ("If the interactive window is run from the wrong directory, it won't work. ")
    print ("Do you need from . ? ")


# To run this code in jupyter, click play. 
# Then create the object in the jupyter terminal. 

class threed_models():
    # This class can fit either the Jorgensen or CMEM model with a Nelder-Mead minimisation routine.  
    # It can also switch between three different cost functions. 

    def __init__(self, filename="S05D05V400B0000-05rad.fits", ppmlr=None, \
                  xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None):
        # This can take arrays of data into the function so real data 
        # can be modelled. Default is it's blank so test data will be used. 
        
        
        if ppmlr is None: # Read in the data from the simulation file if not passed in. 
            self.filename = filename
            ts = process_time()
            print ("Reading ppmlr data:")
            #ppmlr = read_ppmlr.read_ppmlr_cube(filename=self.filename, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
            ppmlr = read_fits_cube.read_fits_cube(filename=self.filename, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
            te = process_time()
            print ("Time = {:.1f}s".format(te-ts))

        else:
            self.filename = ppmlr.filename.split("/")[-1]

        # Extract the x, y, z and eta data from the file object. 
        # self.ppmlr = ppmlr
        self.n = ppmlr.n
        self.x = ppmlr.x_3d
        self.y = ppmlr.y_3d
        self.z = ppmlr.z_3d
        self.eta = ppmlr.eta_3d 

        # Extract any useful solar wind parameters
        self.temp = ppmlr.temp
        self.density = ppmlr.density
        self.vx = ppmlr.vx
        self.vy = ppmlr.vy
        self.vz = ppmlr.vz
        self.bx = ppmlr.bx
        self.by = ppmlr.by
        self.bz = ppmlr.bz
        self.pdyn = ppmlr.dyn_pressure
        self.pmag = ppmlr.mag_pressure
        self.dipole = 0 

        # Get the r, theta and phi coordinates. 
        # Convert to Shue cooords. to calculate the function. 
        print ("Calculating shue coordinates:")
        ts = process_time()
        self.r, self.theta, self.phi = cconv.convert_xyz_to_shue_coords(self.x, self.y, self.z)
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))   
        
        

    def __repr__(self):
        # This is what will be printed about the model if you print the object to the terminal. 
        # It's meant to be better than using __string__(self):
        return f"3D model object. Can do 'jorg' or 'cmem' "


    #COST FUNCTIONS
    ###############
    
    def get_cost_function(self):
        '''This returns the cost function that calculates the misfit/n.
        
        Parameters
        ----------
        self - variable that contains the data.  
        
        Returns
        -------
        Cost Function. 
            - if self.cost_func == "sum squares", it will return the cost function using squared deviations.  
            - elif self.cost_func == "absolute", it will return the cost function using absolute deviations. 
         
        '''

        # To minimise with Nelder-Mead, you need a cost 
        #function, or as Matt J. called it, the misfit 
        #function. See his FitPlasma.py file in WaveHarmonics. 

        if self.cost_func.lower() == "sum squares":
            def cost_func_sum_squared_differences_by_n(params):
                # This cost function is the sum of 
                #the squared differences divided by n. 
                
                # Calculate Model values of eta with 
                #a given set of parameters.
                eta_model = self.get_eta_model(params) 
                    
                # Now get the misfit, (model - observed)**2
                sq_diff = (eta_model - self.eta)**2
                cost = sq_diff.sum()/self.eta.size
                self.cost_per_iteration.append(cost)
                self.param_list.append(params)
                print (cost)
               
                return cost
            return cost_func_sum_squared_differences_by_n
        
        elif self.cost_func.lower() == "absolute":
            def cost_func_sum_absolute_deviations_by_n(params):
                # This cost function is the sum of 
                #the absolute deviations divided by n. 
                # This is the cost function used in 
                #Jorgensen et al. (2019b).

                # Calculate Model values of eta with a 
                #given set of parameters.
                eta_model = self.get_eta_model(params) 
                

                #  Now get the misfit, abs(model - observed)
                abs_diff = abs(eta_model - self.eta)
                cost = abs_diff.sum()/self.eta.size
                self.cost_per_iteration.append(cost)
                self.param_list.append(params)
                print (cost)

                return cost
            return cost_func_sum_absolute_deviations_by_n
            
        elif self.cost_func.lower() == "normalised":
            def cost_func_sum_squares_by_sum_observed(params):
                # This cost function is the sum of the 
                #squared deviations normalised by the data value
                # and the sum of the observed emissivities. 

                # Calculate Model values of eta with 
                #a given set of parameters.
                eta_model = self.get_eta_model(params) 
                
                
                #  Now get the misfit, (model - observed)**2
                sq_diff_norm = (eta_model - self.eta)**2
                cost = sq_diff_norm.sum()/(self.eta**2).sum()
                self.cost_per_iteration.append(cost)
                self.param_list.append(params)
                print (cost)

                return cost
            return cost_func_sum_squares_by_sum_observed
        else:
            raise ValueError("Invalid cost function chosen. Select either 'sum squares', 'absolute' or 'normalised'.") 

    def get_eta_model(self, params):
        '''This function calculates the eta model values for each iteration. This function is intended to be run from the cost function. 
        Parameters
        ----------
        params - tuple of the model parameters for the chosen model. '''
        
        if self.current_model == "jorg": 
            eta_model = self.current_func(self.r, self.theta, self.phi, *params)
        elif self.current_model == "cmem":
            eta_model = self.current_func(self.r, self.theta, self.phi, *self.lin_coeffs, *params)
        else: 
            raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))
        
        return eta_model


        
    #FITTING THE MODEL TO THE DATA
    ##############################
           
    def fit_function_with_nelder_mead(self, model = "jorg", params0 = None, set_param_bounds=False, cost_func="normalised", init_method=1):
        '''This uses a Nelder-Mead minimisation technique to find the best 
        parameters for the chosen model to the data. 
        
        Parameters
        ----------
        model - which model to fit. "jorg" or "cmem" 
        params - (a0, b0, ...) - Tuple containing initial guesses for the model parameter.
            def = None. Will default to values inside the program for each model unless specified. 
              b0 - Initial guess for the intercept parameter b 
        param_bounds - Boolean. If true, it will apply boundaries on the values of the parameters to be fitted. 
        cost - Type of cost function to use. "sum squares" (def), "absolute" or "normalised"
            - Sum Squares will calculate the sum of the squared deviations/n
            - Absolute will calculate the sum of the absolute deviations/n 
            - Normalised will calculate the sum of the (squared deviations) /(n*sum observed)
        init_method - Boolean to use either method 1 or method 2 from the CMEM paper to set the initial model parameters. 
        pickle - boolean to create the pickled file. 

        '''

        
        # Use this tutorial: https://machinelearningmastery.com/how-to-use-nelder-mead-optimization-in-python/#:~:text=The%20Nelder%2DMead%20optimization%20algorithm%20can%20be
        #%20used%20in%20Python,initial%20point%20for%20the%20search.
        # Nelder-Mead does not use gradient methods to find the best fit. 
        # It requires a starting point. 
        # It can be used for multidimensional functions with multiple parameters. 
        # It can be applied to multi-modal functions. 
        # The correct reference is at: https://academic.oup.com/comjnl/article-abstract/7/4/308/354237?login=true

        # First, select the correct model to fit. 
        self.current_model = model.lower() 
        self.current_func = bef.get_model_func(self.current_model)
        self.init_method = init_method

        #GET THE INITIAL PARAMETERS (IF NOT GIVEN)
        ##########################################
        
        # Get Lin model coefficients. 
        if self.current_model == "cmem":
            self.lin_coeffs = bef.get_lin_coeffs(self.dipole, self.pdyn, self.pmag, self.bz)
            self.r0_lin = self.lin_coeffs[-1] 
            
        #Get initial parameters. 
        if params0 is None: 
            if self.current_model == "jorg":
                
                self.params0 = sip.get_init_params(self.current_model, self.init_method, self.bz, self.pdyn, self.density) 
        
            elif self.current_model == "cmem":
                self.params0 = sip.get_init_params(self.current_model, self.init_method, self.bz, self.pdyn, self.density, self.r0_lin) 
        
            else:
                raise ValueError("{} not a valid model. Choose 'cmem' or 'jorg'".format(self.current_model))
        else:
            self.params0 = params0 
        
       


        #SET BOUNDARIES ON PARAMETERS (OPTIONAL - NOT USED)
        ###################################################
        
        #if set_param_bounds: 
            
        #    if self.current_model == "jorg":
                # Set boundaries on some parameters: 
                #MP, BS, alpha parameters.
                # Mp and Bs locations can vary up to 3Re 
                #from initial guess. 
        #        self.param_bounds = ((self.params0[0]-3,self.params0[0]+3), (self.params0[1]-3,self.params0[1]+3), (None,None), (None,None), (None,None), (None,None), (None,None), (0.5,1.5), (0.5,1.5), (0.5,1.5), (0.5,1.5))
        #    elif self.current_model == "cmem":
                # Set boundaries on some parameters. 
                #Default will be not to run this at first. 
        #        self.param_bounds = ((None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None))
        #    else: 
        #        raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))
          
        #else: 
        self.param_bounds=None

        #GET COST FUNCTION AND MINIMISE IT.
        ###################################

        # Get cost calculation function. 
        self.cost_func = cost_func.lower()
        Calculate_cost = self.get_cost_function()
        
        # Set arrays to record info as optimisation progresses. 
        self.cost_per_iteration = []
        self.param_list = []
        

        # The minimize function takes the function and 
        #initial parameters to search.
        # There is an option to add in boundaries for each 
        #parameter as (min, max) pairs.
        # It returns an OptimizeResult object (see scipy docs). 
        print ("Minimising function:")
        ts = process_time() 
        self.result = minimize(Calculate_cost, self.params0, method='nelder-mead', bounds=None)
        te = process_time()
        self.opt_time = te-ts
        print ("Time = {:.1f}s".format(self.opt_time)) 

        # Extract the best fit parameters of a and b, 
        #as well as the final cost and number of iterations. 
        self.params_best_nm = self.result.x
        self.minimum_cost = self.result.fun 
        
        # This is not the number of times it went 
        #through the cost function... 
        self.iterations = self.result.nfev 
        self.cost_per_iteration = np.array(self.cost_per_iteration)

        self.eta_model = self.get_eta_model(self.params_best_nm) 
  

    def write_pickle(self, fname=None, savetag=""):
        '''This will create a pickle file of all the information that would be needed for plotting.
        This is to save an object already created. 
        
        
        '''

        # Name and locate the pickle file. 
        pickle_path = os.environ.get("PICKLE_PATH")
        
        if fname is None:
            fname = self.filename+"_{}_{}_im{}_{}.pkl".format(self.current_model, self.cost_func, self.init_method, savetag)

        # Add all the desired information to a dictionary. 
        pickle_dict = {
                       "cost func":self.cost_func,
                       "min cost":self.minimum_cost,
                       "param list":self.param_list,
                       "cost per it":self.cost_per_iteration,
                       "param bounds":self.param_bounds,
                       "opt time":self.opt_time,
                       "x":self.x,
                       "y":self.y,
                       "z":self.z,
                       "etad":self.eta,
                       "etam":self.eta_model,
                       "params0":self.params0,
                       "params best nm":self.params_best_nm,
                       "filename":self.filename,
                       "model":self.current_model,
                       "r0lin":self.r0_lin if self.current_model == "cmem" else 0,
                       "temp":self.temp,
                       "density":self.density,
                       "vx":self.vx,
                       "vy":self.vy,
                       "vz":self.vz,
                       "bx":self.bx,
                       "by":self.by,
                       "bz":self.bz,
                       "pdyn":self.pdyn,
                       "pmag":self.pmag,
                       "dipole":self.dipole,
                       "init method":self.init_method
                       }
        
        with open(os.path.join(pickle_path, self.current_model+"_optimised", fname), "wb") as f:
            pickle.dump(pickle_dict, f)

        print ("Pickled: {}".format(fname))

