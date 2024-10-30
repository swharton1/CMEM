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
                 xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, model="jorg"):

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
        if params0 is None:
            if self.current_model == "jorg":
                 # Get initial Mp using Shue et al. (1997) formula. Initial Bs is Mp + 3. 
                mp_i = self.get_initial_magnetopause() 

                # Get initial alpha values for Mp. Bs values are Mp + 0.2. 
                alpha_i = self.get_initial_alpha()

                self.params0 = (mp_i,mp_i+3,0.000032, 0.000013, -0.000018, 2.5, -1.6, alpha_i, alpha_i, alpha_i+0.2, alpha_i+0.2)
                # params0 = (8, 11, 0.000032, 0.000013, -0.000018, 2.5, -1.6, 0.6, 0.4, 0.8, 0.8)
            
            elif self.current_model == "cmem":
                self.params0 = (1,12,0.000015, 0.000013, 2, 2.5, -1.6, 1, 3, 4, 0.8, 0.8)
            else:
                raise ValueError("{} not a valid model. 'jorgensen' or 'cmem' only atm.".format(self.current_model))
        else:
            self.params0 = params0
        print (self.params0)

        # Now calculate eta from the model. 
        print ("Calculating eta with model: ")
        ts = process_time()
        
        self.current_func = self.get_model_func()

        self.eta_model = self.current_func(self.r, self.theta, self.phi, *self.params0)
        te = process_time()
        print ("Calculated model eta values: {:.1f}s".format(te-ts))

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
            
#    def convert_xyz_to_shue_coords(self, x, y, z):
#        '''This will convert the x,y,z coordinates to those used in the Shue model 
#         of the magnetopause and bowshock. 

#        Parameters
#        ----------
#        x, y, z - now 3D.  

#        Returns
#        -------
#        r, theta (rad) and phi (rad)
#        '''

#        # r 
#        r = (x**2 + y**2 + z**2)**0.5
        
#        # theta - only calc. where coordinate singularities won't occur. 
#        theta = np.zeros(r.shape)
#        i = np.where(r != 0)
#        theta[i] =  np.arccos(x[i]/r[i])

        # phi - only calc. where coordinate singularities won't occur. 
#        phi = np.zeros(r.shape)
#        j = np.where((y**2 + z**2) != 0)
#        phi[j] = np.arccos(y[j]/((y[j]**2 + z[j]**2)**0.5))
        
#        return r, theta, phi
    
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
   
    def lin_scaled_func(self, theta, phi, dipole=0, pd=2, pm=0.0001, bz=-0.5, p0=1, p1=1, p2=1, p3=1):
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
            p2 scales indentation parameter Q 
            p3 scales d in indentation shape
            '''

        # Get a coefficients first. 
        a = np.array([12.544, -0.194, 0.305, 0.0573, 2.178, 0.0571, -0.999, 16.473, 0.00152, 0.382, 0.0431, -0.00763, -0.210, 0.0405, -4.430, -0.636, -2.600, 0.832, -5.328, 1.103, -0.907, 1.450])

        # Get beta coefficients. 
        beta = np.array([a[6] + a[7]*((np.exp(a[8]*bz) - 1)/(np.exp(a[9]*bz) + 1)), a[10], a[11] + a[12]*dipole, a[13]])

        # Get cn and cs coefficients (equal). 
        c = a[14]*(pd+pm)**a[15]

        # Get d coefficients. 
        dn = p3*(a[16] + a[17]*dipole + a[18]*dipole**2)
        ds = p3*(a[16] - a[17]*dipole + a[18]*dipole**2)
        
        # Get theta-n and theta-s coefficients.
        theta_n = a[19] + a[20]*dipole
        theta_s = a[19] - a[20]*dipole
        
        # Get phi-n and phi-s.
        phi_n = np.arccos((np.cos(theta)*np.cos(theta_n)) + (np.sin(theta)*np.sin(theta_n)*np.cos(phi-(np.pi/2.))))
        phi_s = np.arccos((np.cos(theta)*np.cos(theta_s)) + (np.sin(theta)*np.sin(theta_s)*np.cos(phi-(3*np.pi/2.))))

        # Get f. 
        f = (np.cos(theta/2) + a[5]*np.sin(2*theta)*(1-np.exp(-theta)))**(p1*(beta[0] + beta[1]*np.cos(phi) + beta[2]*np.sin(phi) + beta[3]*(np.sin(phi)**2)))

        # Get r0. 
        self.r0_lin = a[0]*((pd+pm)**a[1])*(1 + a[2]*((np.exp(a[3]*bz) -1 )/(np.exp(a[4]*bz) + 1)))
        print (self.r0_lin)
        # Get Q. 
        Q = p2*c*np.exp(dn*(phi_n**a[21])) + p2*c*np.exp(ds*(phi_s**a[21]))

        # Get r. 
        r = p0*self.r0_lin*f + Q

        return r

    def get_model_func(self):
        '''This will select the correct model function. '''
        if self.current_model == "jorg":
            def jorgensen_func(r, theta, phi, mp, bs, A1, A2, B, alpha, beta, ay_mp, az_mp, ay_bs, az_bs):
               
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
            return jorgensen_func
       
        
        elif self.current_model == "cmem":
            def jorgensen_mod_lin_func(r, theta, phi, p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs):
                '''
                This is the adapted jorgensen-mod model, which will use the lin model to work out 
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
            return jorgensen_mod_lin_func
        
        else:
            raise ValueError("{} not a valid model. 'jorgensen' or 'cmem' only atm.".format(self.current_model))

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

        
