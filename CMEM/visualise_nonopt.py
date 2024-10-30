import numpy as np 
from matplotlib.patches import Wedge, Polygon, Circle
from time import process_time
import matplotlib.pyplot as plt
import os 

# This function will compare ppmlr files with a 
# jorgensen model with stated parameters visually. 
try: 
    #from . import read_ppmlr
    from . import ppmlr_fits
    from . import get_names_and_units as gnau 
    from . import get_meridians as gm 
    from . import coord_conv as cconv 
    from . import boundary_emissivity_functions as bef
    from . import set_initial_params as sip 
    
    
except(ImportError):
    print ("Are you working from the right directory? ")
    print ("If the interactive window is run from the wrong directory, it won't work. ")
    print ("Do you need from . ? ")

class compare_data_model():
    '''This function will read in the PPMLR cube files to get the emissivity 
    data from model runs. It will also create a model with the stated 
    parameters and plot them side by side to visually compare them. 
    '''

    def __init__(self, filename="S05D05V400B0000-05rad.fits", ppmlr=None, params0 = None, \
                 xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, model="jorg", init_method=1):

        # Read in the data from the simulation file. 
        if ppmlr is None: 
            ts = process_time()
            print ("Reading ppmlr data:")
            self.filename=filename
            #ppmlr = read_ppmlr.read_ppmlr_cube(filename=self.filename, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
            ppmlr = ppmlr_fits.read_ppmlr_fits(filename=self.filename, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
            #ppmlr.reshape_to_3D()
            #ppmlr.apply_limit(c)
            te = process_time()
            print ("Completed read in of ppmlr data: {:.1f}s".format(te-ts))
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
        self.density = ppmlr.density
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
        print ("Calculated shue coordinates: {:.1f}s".format(te-ts)) 


        # Now record the parameters you wish to put into the model. 
        self.current_model = model.lower() 
        self.init_method = init_method
        
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
        print (self.params0)

        # Now calculate eta from the model. 
        print ("Calculating eta with model: ")
        ts = process_time()
        
        self.current_func = bef.get_model_func(self.current_model)

        if self.current_model == "jorg": 
            self.eta_model = self.current_func(self.r, self.theta, self.phi, *self.params0)
        elif self.current_model == "cmem":
            self.eta_model = self.current_func(self.r, self.theta, self.phi, *self.lin_coeffs, *self.params0)
        else: 
            raise ValueError("{} not a valid model. 'jorg' or 'cmem' only atm.".format(self.current_model))


        if self.current_model == "cmem":
            self.image_tag = "CMEM"
        else:
            self.image_tag = self.current_model.capitalize()



    def __repr__(self):
        return ("Custom compare_data_model object for the file: {}".format(self.filename))
    
    def modify_parameters_jorg(self, mp=None, bs=None, A1=None, A2=None): 
        '''This function will manually update the parameters mp, bs, A1 and A2. Use for CMEM model. '''
        print (mp, bs, A1, A2) 
        
        #Update parameters. 
        new_params = (mp, bs, A1, A2, *self.params0[4:])
        self.params0 = new_params  
        
    def modify_parameters_cmem(self, p0=None, bs=None, A1=None, A2=None): 
        '''This function will manually update the parameters p0, bs, A1 and A2. Use for CMEM model. '''
        print (p0, bs, A1, A2) 
        
        #Update parameters. 
        new_params = (p0, bs, A1, A2, *self.params0[4:])
        self.params0 = new_params 
            
    #PLOTTING FUNCTIONS
    ###################
          
    def plot_planes(self, cmap='hot', vmin=-8, vmax=-4, levels=100, save=False, savetag=""):
        '''This will just plot the x-z and x-y planes through the model (recommended way).
        
        Parameters
        ----------
        cmap - matplotlib colourmap.
        vmin - minimum log value of eta to show on the contour map. def = -8
        vmax - maximum log value of eta to show on the contour map. def = -3
        levels - number of levels on the contour map. 
        '''
 
        #Get meridian data for etad. 
        xp_y, yp_y, zp_y, etad_y, xp_z, yp_z, zp_z, etad_z = gm.calculate_meridian_planes(self.x, self.y, self.z, self.eta)
        
        #Get meridian data for etam. 
        xp_y, yp_y, zp_y, etam_y, xp_z, yp_z, zp_z, etam_z = gm.calculate_meridian_planes(self.x, self.y, self.z, self.eta_model)
        
        
        # Calculate log10 eta values. If eta = 0, set log(eta) = -12  
        letad_y = np.zeros(etad_y.shape)+vmin
        i = np.where(etad_y != 0)
        letad_y[i] = np.log10(etad_y[i])
        j = np.where(letad_y < vmin)
        letad_y[j] = vmin 

        letam_y = np.zeros(etam_y.shape)+vmin
        i = np.where(etam_y != 0)
        letam_y[i] = np.log10(etam_y[i])
        j = np.where(letam_y < vmin)
        letam_y[j] = vmin 

        letad_z = np.zeros(etad_z.shape)+vmin
        i = np.where(etad_z != 0)
        letad_z[i] = np.log10(etad_z[i])
        j = np.where(letad_z < vmin)
        letad_z[j] = vmin 

        letam_z = np.zeros(etam_z.shape)+vmin
        i = np.where(etam_z != 0)
        letam_z[i] = np.log10(etam_z[i])
        j = np.where(letam_z < vmin)
        letam_z[j] = vmin 
        
        # Create a filename label so you know which file you plotted. 
        file_label = self.filename.split("/")[-1]

        # Now you can make the contour plot. 
        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(bottom=0.20, hspace=0.4, wspace=0.2)

        # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)

        # etad_y
        ax1 = fig.add_subplot(221)
        cont1 = ax1.contourf(xp_y, zp_y, letad_y, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        ax1.set_xlabel('X [RE]')
        ax1.set_ylabel('Z [RE]')
        ax1.set_title("n = {:.2f} cm".format(self.density)+r"$^{-3}$"+"\nXZ Plane")
        ax1.set_aspect("equal")
        self.make_earth(ax1, rotation=-90)

        self.cont1 = cont1

        # Colourbars 
        cbar = plt.colorbar(cont1, ax=ax1)
        # cbar.set_label(r"eV cm$^{-3}$ s$^{-1}$")
        level_min = int(np.ceil(cont1.levels.min()))
        level_max = int(np.floor(cont1.levels.max()))
        cticks = np.arange(level_min, level_max+1)
        cbar.set_ticks(cticks)
        cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])

        # etam_y
        ax2 = fig.add_subplot(222)
        cont2 = ax2.contourf(xp_y, zp_y, letam_y, cmap='hot', levels=cont1.levels, vmin=vmin, vmax=vmax)
        ax2.set_xlabel('X [RE]')
        ax2.set_ylabel('Z [RE]')
        ax2.set_title("{}\nXZ Plane".format(self.image_tag))
        ax2.set_aspect("equal")
        self.make_earth(ax2, rotation=-90)

        # Colourbars 
        cbar = plt.colorbar(cont2, ax=ax2)
        cbar.set_label(r"eV cm$^{-3}$ s$^{-1}$")
        level_min = int(np.ceil(cont2.levels.min()))
        level_max = int(np.floor(cont2.levels.max()))
        cticks = np.arange(level_min, level_max+1)
        cbar.set_ticks(cticks)
        cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])

        # etad_z
        ax3 = fig.add_subplot(223)
        cont3 = ax3.contourf(xp_z, yp_z, letad_z, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        ax3.set_xlabel('X [RE]')
        ax3.set_ylabel('Y [RE]')
        ax3.set_title("XY Plane")
        ax3.set_aspect("equal")
        self.make_earth(ax3, rotation=-90)

        # Colourbars 
        cbar = plt.colorbar(cont3, ax=ax3)
        # cbar.set_label(r"eV cm$^{-3}$ s$^{-1}$")
        level_min = int(np.ceil(cont3.levels.min()))
        level_max = int(np.floor(cont3.levels.max()))
        cticks = np.arange(level_min, level_max+1)
        cbar.set_ticks(cticks)
        cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])

        # etam_z
        ax4 = fig.add_subplot(224)
        cont4 = ax4.contourf(xp_z, yp_z, letam_z, cmap='hot', levels=cont3.levels, vmin=vmin, vmax=vmax)
        ax4.set_xlabel('X [RE]')
        ax4.set_ylabel('Y [RE]')
        ax4.set_title("XY Plane")
        ax4.set_aspect("equal")
        self.make_earth(ax4, rotation=-90)

        # Colourbars 
        cbar = plt.colorbar(cont4, ax=ax4)
        cbar.set_label(r"eV cm$^{-3}$ s$^{-1}$")
        level_min = int(np.ceil(cont4.levels.min()))
        level_max = int(np.floor(cont4.levels.max()))
        cticks = np.arange(level_min, level_max+1)
        cbar.set_ticks(cticks)
        cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])
        
        # Add a label to show the model parameters. 
        label = ""
        info = gnau.get_parameter_info(model=self.current_model)
        parameter_names = [info[i][0] for i in info.keys()]
        parameter_units = [info[i][1] for i in info.keys()]
        for p,pval in enumerate(self.params0):
                pv = pval 
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pv,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"

        fig.text(0.5, 0.02, label, ha='center')

        self.fig = fig 

        if save: 
            fig.savefig(os.environ.get("PLOT_PATH")+"{}/{}_data_{}_model_planes_nonopt{}.png".format(self.current_model, self.filename, self.current_model, savetag))
           

    def plot_earth_sun_line(self, save=False, savetag=""):
        '''This will plot the emissivity along a line close to the sun-earth line. 
        This will make it easier to see how the function compares to the simulation. 
        '''

        #Get Earth_sun line data for emissivity data. 
        xp, yp, zp, etad = gm.calculate_sunearth_line(self.x, self.y, self.z, self.eta)
        
        #Get Earth_sun line data for emissivity model. 
        xp, yp, zp, etam = gm.calculate_sunearth_line(self.x, self.y, self.z, self.eta_model)


        # Separate the model line into three colours for the different model sections. 
        if self.current_model == "jorg":
            i_msphere = np.where(xp <= self.params0[0])
            i_msheath = np.where((xp > self.params0[0]) & (xp <= self.params0[1]))
            i_bow = np.where(xp > self.params0[1])
        elif self.current_model == "cmem":
            i_msphere = np.where(xp <= self.params0[0]*self.r0_lin)
            i_msheath = np.where((xp > self.params0[0]*self.r0_lin) & (xp <= self.params0[1]))
            i_bow = np.where(xp > self.params0[1])
            
        # Now you can make the contour plot. 
        fig = plt.figure(figsize=(8,6))
        fig.subplots_adjust(bottom=0.20, top=0.80)
        ax = fig.add_subplot(111)

        ax.plot(xp, etad, 'k', label="PPMLR")
        ax.plot(xp[i_msphere], etam[i_msphere], 'b', label="Model - Magnetosphere")
        ax.plot(xp[i_msheath], etam[i_msheath], 'r', label="Model - Magnetosheath")
        ax.plot(xp[i_bow], etam[i_bow], 'g', label="Model - Solar Wind")
        
        ax.set_xlabel('X [RE]')
        ax.set_ylabel(r"eV cm$^{-3}$ s$^{-1}$")
        ax.legend(loc='best')
        ax.set_title("Simulation Data vs {} Model - Sun-Earth Line\nn = {:.2f} cm".format(self.image_tag, self.density)+r"$^{-3}$"+"\nInitial Parameters")
        ax.set_xlim(0,25)

    # Add a label to show the model parameters. 
        label = ""
        info = gnau.get_parameter_info(model=self.current_model)
        parameter_names = [info[i][0] for i in info.keys()]
        parameter_units = [info[i][1] for i in info.keys()]
        for p,pval in enumerate(self.params0):
                pv = pval 
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pv,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"
        

        fig.text(0.5, 0.02, label, ha='center')
        self.fig_sunearth = fig 

        if save: 
            print ("{}/{}_data_{}_model_sunearth_nonopt{}.png".format(self.current_model, self.filename, self.current_model, savetag))
            fig.savefig(os.environ.get("PLOT_PATH")+"{}/{}_data_{}_model_sunearth_nonopt{}.png".format(self.current_model, self.filename, self.current_model, savetag))
           




    def make_earth(self, ax, rotation=0):
        '''This will add a little plot of the Earth on top for reference. '''

        # Add white circle first. 
        r=1
        circle = Circle((0,0), r, facecolor='w', edgecolor='navy')
        ax.add_patch(circle)

        # Add nightside. 
        theta2 = np.arange(181)-180+rotation
        xval2 = np.append(r*np.cos(theta2*(np.pi/180)),0)
        yval2 = np.append(r*np.sin(theta2*(np.pi/180)),0)
        verts2 = [[xval2[i],yval2[i]] for i in range(len(xval2))]
        
        polygon2 = Polygon(verts2, closed=True, edgecolor='navy', facecolor='navy', alpha=1) 
        ax.add_patch(polygon2)
    
    def sig_figs(self, x: float, precision: int):
        """
        Rounds a number to number of significant figures
        Parameters:
        - x - the number to be rounded
        - precision (integer) - the number of significant figures
        Returns:
        - float
        """

        x = float(x)
        precision = int(precision)

        return np.round(x, -int(np.floor(np.log10(abs(x)))) + (precision - 1))

        
