import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import string 
from time import process_time
from matplotlib.patches import Wedge, Polygon, Circle
import pickle
import os

#This is the original version of the fitting code that has all the 
#functions stored internally instead of reading them from 
#boundary_emissivity_functions.py. I kept it in case I screwed up 
#the code transfer! 

try: 
    #from . import read_ppmlr
    from . import ppmlr_fits
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
            ppmlr = ppmlr_fits.read_ppmlr_fits(filename=self.filename, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
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
        self.r, self.theta, self.phi = self.convert_xyz_to_shue_coords(self.x, self.y, self.z)
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))   
        
        

    def __repr__(self):
        # This is what will be printed about the model if you print the object to the terminal. 
        # It's meant to be better than using __string__(self):
        return f"3D model object. Can do 'jorg' or 'cmem' "

    
    def convert_xyz_to_shue_coords(self, x, y, z):
        '''This will convert the x,y,z coordinates to those used in the Shue model 
         of the magnetopause and bowshock. 

        Parameters
        ----------
        x, y, z - now 3D.  

        Returns
        -------
        r, theta (rad) and phi (rad)
        '''

        # r 
        r = (x**2 + y**2 + z**2)**0.5
        
        # theta - only calc. where coordinate singularities won't occur. 
        theta = np.zeros(r.shape)
        i = np.where(r != 0)
        theta[i] =  np.arccos(x[i]/r[i])

        # phi - only calc. where coordinate singularities won't occur. 
        phi = np.zeros(r.shape)
        j = np.where((y**2 + z**2) != 0)
        phi[j] = np.arccos(y[j]/((y[j]**2 + z[j]**2)**0.5))
        
        return r, theta, phi
        
    def convert_shue_to_xyz_coords(self, r, theta, phi):
        '''This will convert the Shue coordinates back to xyz coordinates. 
        
        Parameters
        ----------
        r, theta (rad), phi (rad)
        
        Returns
        -------
        x,y,z
        '''

        x = r*np.cos(theta)
        y = r*np.sin(theta)*np.cos(phi)
        z = r*np.sin(theta)*np.sin(phi)

        return x,y,z 

    def shue_func(self, theta, phi, r0, ay, az):
        '''This is the 3D Shue model defined in Jorgensen et al. (2019)
        
        Parameters
        ----------
        theta (rad) and phi (rad)
        r0 - subsolar magnetopause distance
        ay - alpha y parameter
        az - alpha z parameter 

        Returns
        -------
        r - radial distance at the angles theta and phi 
        '''

        ry = r0*((2/(1+np.cos(theta)))**ay)
        rz = r0*((2/(1+np.cos(theta)))**az)

        r = (ry*rz)/(((rz*np.cos(phi))**2 + (ry*np.sin(phi))**2)**0.5)

        return r 

    def lin_scaled_func(self, theta, phi, dipole=0, pd=2, pm=0.01, bz=-0.5, p0=1, p1=1, p2=1, p3=1):
        '''This function will work out r using the lin model. 
        
        Parameters
        ----------
        theta (rad) - Shue coords.
        phi (rad) - Shue coords. 
        dipole - dipole tilt angle (rad)
        pd - dynamic pressure in nPa
        pm - magnetic pressure in nPa 
        bz - IMF bz component in nT 
        p - parameter scaling factors. 
            p0 scales r0
            p1 scales flaring parameter beta 
            p2 scales indentation parameter Q (cusp depth) 
            p3 scales d in indentation shape (cusp shape/width)
            '''

        # Get coefficients if for some reason, they have not already been calculated. 
        if self.r0_lin is None: 
            self.get_lin_coeffs(dipole, pd, pm, bz)
        
        # Get phi-n and phi-s.
        phi_n = np.arccos((np.cos(theta)*np.cos(self.theta_n)) + (np.sin(theta)*np.sin(self.theta_n)*np.cos(phi-(np.pi/2.))))
        phi_s = np.arccos((np.cos(theta)*np.cos(self.theta_s)) + (np.sin(theta)*np.sin(self.theta_s)*np.cos(phi-(3*np.pi/2.))))

        # Get f. 
        f = (np.cos(theta/2) + self.a[5]*np.sin(2*theta)*(1-np.exp(-theta)))**(p1*(self.beta_c[0] + self.beta_c[1]*np.cos(phi) + self.beta_c[2]*np.sin(phi) + self.beta_c[3]*(np.sin(phi)**2)))

        # Get Q. 
        Q = p2*self.c*np.exp(self.dn*(phi_n**self.a[21])) + p2*self.c*np.exp(self.ds*(phi_s**self.a[21]))

        # Get r. 
        r = p0*self.r0_lin*f + Q

        return r 


    def get_model_func(self):
        '''This will select the correct function for the desired model. '''
        
        if self.current_model == "jorg":
            def jorg_func(r, theta, phi, mp, bs, A1, A2, B, alpha, beta, ay_mp, az_mp, ay_bs, az_bs):
               
                '''This is the model from the Jorgensen paper. 
        
                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                mp - subsolar magnetopause distance parameter
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                ay_mp - ay magnetopause flaring parameter
                az_mp - az magnetopause flaring parameter
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                '''

                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = self.shue_func(theta, phi, mp, ay_mp, az_mp)
                rbs = self.shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)

                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (A1 + B*((np.sin(theta[r2]))**8))*((r[r2]/10)**(-alpha-(beta*(np.sin(theta[r2]))**2)))
                eta[r3] = A2*((r[r3]/10)**(-3))
        
                return eta
            return jorg_func
             
        elif self.current_model == "cmem":
            def cmem_func(r, theta, phi, p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs):
                '''
                This is the CMEM model, which will use the lin model to work out 
                the magnetopause location instead of the shue model. 

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                p0 - scaling factor on the subsolar magnetopause parameter 
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                p1 - scaling factor on magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                '''
            
                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = self.lin_scaled_func(theta, phi, 0, 2, 0.01, -0.5, p0, p1, p2, p3)
                rbs = self.shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)

                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = A1*(np.exp(-B*(theta[r2]/2.)**4))*((r[r2]/10)**(-alpha-(beta*(np.sin(theta[r2]))**2)))
                eta[r3] = A2*((r[r3]/10)**(-3))
                
                return eta
            return cmem_func
        
        else:
            raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))


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

        # To minimise with Nelder-Mead, you need a cost function, or as Matt J. called it,
        # the misfit function. See his FitPlasma.py file in WaveHarmonics. 

        if self.cost_func.lower() == "sum squares":
            def cost_func_sum_squared_differences_by_n(params):
                # This cost function is the sum of the squared differences divided by n. 
                
                # Calculate Model values of eta with a given set of parameters.
                eta_model = self.current_func(self.r, self.theta, self.phi, *params)
                
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
                # This cost function is the sum of the absolute deviations divided by n. 
                # This is the cost function used in Jorgensen et al. (2019b).

                # Calculate Model values of eta with a given set of parameters.
                eta_model = self.current_func(self.r, self.theta, self.phi, *params)
                

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
                # This cost function is the sum of the squared deviations normalised by the data value
                # and the sum of the observed emissivities. 

                # Calculate Model values of eta with a given set of parameters.
                eta_model = self.current_func(self.r, self.theta, self.phi, *params)
                
                
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

    #METHOD 1 INITIALISING FUNCTIONS 

    def get_initial_magnetopause(self):
        '''This uses equation 12 in Shue et al. (1997) to estimate the initial 
        subsolar magnetopause position from Bz and Dp, which are both in the ppmlr object. '''

        if self.bz >= 0: 
            return (11.4 + 0.013*self.bz)*(self.pdyn**(-1.0/6.6))
        else:
            return (11.4 + 0.14*self.bz)*(self.pdyn**(-1.0/6.6))

    def get_initial_alpha(self):
        '''This uses equation 13 in Shue et al. (1997) to estimate the initial value of alpha
        in the Shue magnetopause model. Assumes all initial alpha values will be equal to start with.'''

        return (0.58 - 0.010*self.bz)*(1 + 0.010*self.pdyn)
    
    def get_lin_coeffs(self, dipole, pd, pm, bz):
        '''This gets the value of r0 in the Lin et al. (2010) model, which is a constant value 
        that depends on solar wind parameters. All of these functions are independent of beta and gamma. 
        
        Parameters
        ----------
        dipole
        pd
        pm
        bz
        
        Returns
        -------
        All coefficients are attached to self. 
        '''

        # Get a coefficients first. 
        a = np.array([12.544, -0.194, 0.305, 0.0573, 2.178, 0.0571, -0.999, 16.473, 0.00152, 0.382, 0.0431, -0.00763, -0.210, 0.0405, -4.430, -0.636, -2.600, 0.832, -5.328, 1.103, -0.907, 1.450])
        self.a = a

        # Get beta coefficients - renamed delta. 
        self.beta_c = np.array([a[6] + a[7]*((np.exp(a[8]*bz) - 1)/(np.exp(a[9]*bz) + 1)), a[10], a[11] + a[12]*dipole, a[13]])
         
        # Get cn and cs coefficients (equal). 
        self.c = a[14]*(pd+pm)**a[15]

        # Get d coefficients. 
        self.dn = (a[16] + a[17]*dipole + a[18]*dipole**2)
        self.ds = (a[16] - a[17]*dipole + a[18]*dipole**2)
        
        # Get theta-n and theta-s coefficients.
        self.theta_n = a[19] + a[20]*dipole
        self.theta_s = a[19] - a[20]*dipole

        # Get the unscaled subsolar magnetopause radius. 
        self.r0_lin = 12.544*((pd+pm)**-0.194)*(1 + 0.305*((np.exp(0.0573*bz) -1 )/(np.exp(2.178*bz) + 1)))


    #METHOD 2 INITIALISING FUNCTIONS 
    
    def get_initial_mp_method2(self, density):
        '''Gets mp for method 2'''
    
        return -0.10*density + 10.28
    
    def get_initial_bs_method2(self, density):
        '''Gets bs for method 2'''
        
        return -0.12*density + 13.24 

    def get_initial_A1_method2(self, density):
        '''This function estimates the initial value of the parameter A1 for the Jorgensen model. '''

        return 0.0000027*density - 0.0000063
    
    def get_initial_A2_method2(self, density):
        '''This function estimates the initial value of the parameter A2. '''

        return 0.0000009*density - 0.0000010

    def get_initial_p0_method2(self, density):
        '''Gets p0 for CMEM for method 2'''
        
        return  0.0022*density + 0.7753
        
        
           
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

        
        # Use this tutorial: https://machinelearningmastery.com/how-to-use-nelder-mead-optimization-in-python/#:~:text=The%20Nelder%2DMead%20optimization%20algorithm%20can%20be%20used%20in%20Python,initial%20point%20for%20the%20search.
        # Nelder-Mead does not use gradient methods to find the best fit. 
        # It requires a starting point. 
        # It can be used for multidimensional functions with multiple parameters. 
        # It can be applied to multi-modal functions. 
        # The correct reference is at: https://academic.oup.com/comjnl/article-abstract/7/4/308/354237?login=true

        # First, select the correct model to fit. 
        self.current_model = model.lower() 
        self.current_func = self.get_model_func()
        self.init_method = init_method

        # Sort out the initial parameters if not specified. 
        if params0 == None:
            
            if self.current_model == "jorg":
                # mp, bs, A1, A2, B, alpha, beta, ay_mp, az_mp, ay_bs, az_bs. 
                
                if self.init_method == 1: 
                    # These values are rounded values taken from a model run in Jorgensen et al, except those calculate using the functions below. 
                
                    # Get initial Mp using Shue et al. (1997) formula. Initial Bs is Mp + 3. 
                    mp_i = self.get_initial_magnetopause() 

                    # Get initial alpha values for Mp. Bs values are Mp + 0.2. 
                    alpha_i = self.get_initial_alpha()

                    self.params0 = (mp_i,mp_i+3, 0.000032, 0.000013, -0.000018, 2.5, -1.6, alpha_i, alpha_i, alpha_i+0.2, alpha_i+0.2)
                    
                elif self.init_method == 2: 
                    mp = self.get_initial_mp_method2(self.density)
                    bs = self.get_initial_bs_method2(self.density)
                    A1 = self.get_initial_A1_method2(self.density)
                    A2 = self.get_initial_A2_method2(self.density)
                
                    # Get initial alpha values for Mp. Bs values are Mp + 0.2. 
                    alpha_i = self.get_initial_alpha()

                    self.params0 = (mp, bs, A1, A2, -0.000018, 2.5, -1.6, alpha_i, alpha_i, alpha_i+0.2, alpha_i+0.2)
                # params0 = (8, 11, 0.000032, 0.000013, -0.000018, 2.5, -1.6, 0.6, 0.4, 0.8, 0.8)
                
            elif self.current_model == "cmem":
                # p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs. 
                
                # Get Lin model coefficients. 
                self.get_lin_coeffs(self.dipole, self.pdyn, self.pmag, self.bz)
                
                
                if self.init_method == 1: 
                    # A1, A2, B, alpha and beta are rounded values from Jorgensen et al. (2019).
                    # p values to scale magnetopause are from inspection of an example. (1,1,3,4) 
                    
                    
                    
                    # Initial Bs is Mp + 3. 
                    bs_i = self.r0_lin + 3

                    # Get initial alpha values for Mp. Bs values are Mp + 0.2. 
                    bs_alpha_i = self.get_initial_alpha() + 0.2

                    
                    self.params0 = (1, bs_i, 0.000015, 0.000013, 2, 2.5, -1.6, 1, 3, 4, bs_alpha_i, bs_alpha_i)
                
                elif self.init_method == 2: 
                    p0 = self.get_initial_p0_method2(self.density)
                    bs = self.get_initial_bs_method2(self.density)
                    A1 = self.get_initial_A1_method2(self.density)
                    A2 = self.get_initial_A2_method2(self.density)
                    
                    # Get initial alpha values for Mp. Bs values are Mp + 0.2. 
                    bs_alpha_i = self.get_initial_alpha() + 0.2
                    
                    self.params0 = (p0, bs, A1, A2, 2, 2.5, -1.6, 1, 3, 4, bs_alpha_i, bs_alpha_i)
                    
            else: raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))
        
        print ('Initial parameters are: ', self.params0)




        # Set boundaries on parameters. 
        if set_param_bounds: 
            
            if self.current_model == "jorg":
                # Set boundaries on some parameters: MP, BS, alpha parameters.
                # Mp and Bs locations can vary up to 3Re from initial guess. 
                self.param_bounds = ((self.params0[0]-3,self.params0[0]+3), (self.params0[1]-3,self.params0[1]+3), (None,None), (None,None), (None,None), (None,None), (None,None), (0.5,1.5), (0.5,1.5), (0.5,1.5), (0.5,1.5))
            elif self.current_model == "cmem":
                # Set boundaries on some parameters. Default will be not to run this at first. 
                self.param_bounds = ((None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None))
            else: 
                raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))
          
        else: 
            self.param_bounds=None

        # Get cost calculation function. 
        self.cost_func = cost_func.lower()
        Calculate_cost = self.get_cost_function()
        
        # Set arrays to record info as optimisation progresses. 
        self.cost_per_iteration = []
        self.param_list = []
        

        # The minimize function takes the function and initial parameters to search.
        # There is an option to add in boundaries for each parameter as (min, max) pairs.
        # It returns an OptimizeResult object (see scipy docs). 
        print ("Minimising function:")
        ts = process_time() 
        self.result = minimize(Calculate_cost, self.params0, method='nelder-mead', bounds=self.param_bounds)
        te = process_time()
        self.opt_time = te-ts
        print ("Time = {:.1f}s".format(self.opt_time)) 

        # Extract the best fit parameters of a and b, as well as the final cost and number of iterations. 
        self.params_best_nm = self.result.x
        self.minimum_cost = self.result.fun 
        
        # This is not the number of times it went through the cost function... 
        self.iterations = self.result.nfev 
        self.cost_per_iteration = np.array(self.cost_per_iteration)

        self.eta_model = self.current_func(self.x, self.y, self.z, *self.params_best_nm)
  

    def write_pickle(self, savetag=""):
        '''This will create a pickle file of all the information that would be needed for plotting.
        This is to save an object already created. 
        
        
        '''

        # Name and locate the pickle file. 
        pickle_path = os.environ.get("PICKLE_PATH")
        
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
                       "bx":self.by,
                       "bx":self.bz,
                       "pdyn":self.pdyn,
                       "pmag":self.pmag,
                       "dipole":self.dipole,
                       "init method":self.init_method
                       }
        
        with open(os.path.join(pickle_path, self.current_model+"_optimised", fname), "wb") as f:
            pickle.dump(pickle_dict, f)

        print ("Pickled: {}".format(fname))

