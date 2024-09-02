import numpy as np 
import matplotlib.pyplot as plt 
import os 
import pickle
from matplotlib.patches import Wedge, Polygon, Circle
import string 
from . import get_names_and_units as gnau 
from . import get_meridians as gm 

# Create a separate object to analyse the results of the fitting procedure. 
# This means the optimisation doesn't need running every time you change a plot! 

class analyse_model():
    def __init__(self, filename="S05D05V400B0000-05rad.dat_jorg_normalised.pkl"):
        '''This takes in the threed model object you've created from a pickle file. '''
        
        self.pickle_path = os.environ.get("PICKLE_PATH")
        self.plot_path = os.environ.get("PLOT_PATH")
        self.filename=filename.split("_")[0]
        self.current_model = filename.split("_")[1]
        if self.current_model == "cmem":
        	self.image_tag = "CMEM"
        else:
        	self.image_tag = self.current_model.capitalize()

        # Read in the pickled optimised model. 
        self.model = self.read_pickle(filename, self.current_model)
        

    def __repr__(self):
        return f"analyse model object."
    
    def read_pickle(self, filename, current_model):
        '''This will read a single pickle file. '''

        with open(os.path.join(self.pickle_path, current_model+"_optimised", filename), 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def plot_change_in_parameters(self, save=False, savetag=""):
        '''This will plot how the parameters changed over the course of the
        optimisation procedure. '''

        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(hspace=0.4, wspace=0.4, top=0.85)
        fig.text(0.5, 0.92, "Parameter variation with optimisation\n{}\nOptimisation Time = {:.1f}s\nModel = {}".format(self.filename, self.model['opt time'], self.image_tag.capitalize()), ha="center")
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
        ax2b = ax2.twinx()
        ax3 = fig.add_subplot(323)
        ax4 = fig.add_subplot(324)
        ax5 = fig.add_subplot(325)
        ax6 = fig.add_subplot(326)
        
        if self.model['model'] == "jorg":
            ax1.text(0.05, 1.05, self.model["parameter list"][0], c="r", transform=ax1.transAxes, fontsize=12)
            ax1.text(0.25, 1.05, self.model["parameter list"][1], c="b", transform=ax1.transAxes, fontsize=12)
            ax2.text(0.05, 1.05, self.model["parameter list"][2], c="r", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.25, 1.05, self.model["parameter list"][3], c="b", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.45, 1.05, self.model["parameter list"][4], c="g", transform=ax2.transAxes, fontsize=12)
            ax3.text(0.05, 1.05, self.model["parameter list"][5], c="r", transform=ax3.transAxes, fontsize=12)
            ax3.text(0.25, 1.05, self.model["parameter list"][6], c="b", transform=ax3.transAxes, fontsize=12)
            ax4.text(0.05, 1.05, self.model["parameter list"][7], c="r", transform=ax4.transAxes, fontsize=12)
            ax4.text(0.25, 1.05, self.model["parameter list"][8], c="b", transform=ax4.transAxes, fontsize=12)
            ax5.text(0.05, 1.05, self.model["parameter list"][9], c="r", transform=ax5.transAxes, fontsize=12)
            ax5.text(0.25, 1.05, self.model["parameter list"][10], c="b", transform=ax5.transAxes, fontsize=12)
           
            
        elif self.model['model'] == "cmem":
            ax1.text(0.05, 1.05, self.model["parameter list"][0], c="r", transform=ax1.transAxes, fontsize=12)
            ax1.text(0.25, 1.05, self.model["parameter list"][1], c="b", transform=ax1.transAxes, fontsize=12)
            ax2.text(0.05, 1.05, self.model["parameter list"][2], c="r", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.25, 1.05, self.model["parameter list"][3], c="b", transform=ax2.transAxes, fontsize=12)
            ax2.text(0.85, 1.05, self.model["parameter list"][4], c="g", transform=ax2.transAxes, fontsize=12)
            ax3.text(0.05, 1.05, self.model["parameter list"][5], c="r", transform=ax3.transAxes, fontsize=12)
            ax3.text(0.25, 1.05, self.model["parameter list"][6], c="b", transform=ax3.transAxes, fontsize=12)
            ax4.text(0.05, 1.05, self.model["parameter list"][7], c="r", transform=ax4.transAxes, fontsize=12)
            ax4.text(0.25, 1.05, self.model["parameter list"][8], c="b", transform=ax4.transAxes, fontsize=12)
            ax4.text(0.45, 1.05, self.model["parameter list"][9], c="g", transform=ax4.transAxes, fontsize=12)
            ax5.text(0.05, 1.05, self.model["parameter list"][10], c="r", transform=ax5.transAxes, fontsize=12)
            ax5.text(0.25, 1.05, self.model["parameter list"][11], c="b", transform=ax5.transAxes, fontsize=12)
        
        # Sort cost axis and values. 
        if self.model['cost func'] == "sum squares":
            ax6.text(0.5, 1.05, str(self.model['cost func'])+" : Min Cost = "+str(self.sig_figs(self.model['min cost']*1e11,3))+r"$x10^{-11}$", c="k", transform=ax6.transAxes, ha="center")
            cpi = self.model['cost per it']*1e11
            ax6.set_ylabel(r"$x10^{-11} (eV cm^{-3} s^{-1})^2$", fontsize=10)
        elif self.model['cost func'] =="absolute":
            ax6.text(0.5, 1.05, str(self.model['cost func'])+" : Min Cost = "+str(self.sig_figs(self.model['min cost']*1e7,3))+r"$x10^{-7}$", c="k", transform=ax6.transAxes, ha="center")
            cpi = self.model['cost per it']*1e7
            ax6.set_ylabel(r"$x10^{-7} eV cm^{-3} s^{-1}$", fontsize=10)
        elif self.model['cost func'] =="normalised":
            ax6.text(0.5, 1.05, str(self.model['cost func'])+" : Min Cost = "+str(self.sig_figs(self.model['min cost'],3)), c="k", transform=ax6.transAxes, ha="center")
            cpi = self.model['cost per it']
            ax6.set_ylabel("Cost", fontsize=10)

        ax1.set_xlabel("Iterations", fontsize=10)
        ax2.set_xlabel("Iterations", fontsize=10)
        ax3.set_xlabel("Iterations", fontsize=10)
        ax4.set_xlabel("Iterations", fontsize=10)
        ax5.set_xlabel("Iterations", fontsize=10)
        ax6.set_xlabel("Iterations", fontsize=10)

        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax2.get_xticklabels() + ax2.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax3.get_xticklabels() + ax3.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax4.get_xticklabels() + ax4.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax5.get_xticklabels() + ax5.get_yticklabels()): 
            label.set_fontsize(8)
        for label in (ax6.get_xticklabels() + ax6.get_yticklabels()): 
            label.set_fontsize(8)

        iteration = np.arange(len(self.model['param list']))
        param_list_t = np.array(self.model['param list']).transpose()
       
        
        if self.model['model'] == "jorg":
            ax1.plot(iteration, param_list_t[0], "r")
            ax1.plot(iteration, param_list_t[1], "b")
            ax2.plot(iteration, param_list_t[2]*100000, "r")
            ax2.plot(iteration, param_list_t[3]*100000, "b")
            ax2.plot(iteration, param_list_t[4]*100000, "g")
            ax3.plot(iteration, param_list_t[5], "r")
            ax3.plot(iteration, param_list_t[6], "b")
            ax4.plot(iteration, param_list_t[7], "r")
            ax4.plot(iteration, param_list_t[8], "b")
            ax5.plot(iteration, param_list_t[9], "r")
            ax5.plot(iteration, param_list_t[10], "b")
            ax6.plot(iteration, cpi, "k", label="Cost")
        elif self.model['model'] == "cmem":
            ax1.plot(iteration, param_list_t[0]*self.model['r0_lin'], "r")
            ax1.plot(iteration, param_list_t[1], "b")
            ax2.plot(iteration, param_list_t[2]*100000, "r")
            ax2.plot(iteration, param_list_t[3]*100000, "b")
            ax2b.plot(iteration, param_list_t[4], "g")
            ax3.plot(iteration, param_list_t[5], "r")
            ax3.plot(iteration, param_list_t[6], "b")
            ax4.plot(iteration, param_list_t[7], "r")
            ax4.plot(iteration, param_list_t[8], "b")
            ax4.plot(iteration, param_list_t[9], "g")
            ax5.plot(iteration, param_list_t[10], "r")
            ax5.plot(iteration, param_list_t[11], "b")
            ax6.plot(iteration, cpi, "k", label="Cost")

        # If boundaries were applied to parameters, plot them on. NOT ADAPTED FOR JORGENSEN LIN. 
        if self.model['param bounds'] is not None: 
            pbounds = self.model['param bounds'] 
            ax1.plot([iteration[0], iteration[-1]], [pbounds[0][0], pbounds[0][0]], 'r--')
            ax1.plot([iteration[0], iteration[-1]], [pbounds[0][1], pbounds[0][1]], 'r--')
            ax1.plot([iteration[0], iteration[-1]], [pbounds[1][0], pbounds[1][0]], 'b--')
            ax1.plot([iteration[0], iteration[-1]], [pbounds[1][1], pbounds[1][1]], 'b--')
            ax4.plot([iteration[0], iteration[-1]], [pbounds[7][0], pbounds[7][0]], 'r--')
            ax4.plot([iteration[0], iteration[-1]], [pbounds[7][1], pbounds[7][1]], 'r--')
            ax4.plot([iteration[0], iteration[-1]], [pbounds[8][0], pbounds[8][0]], 'b--')
            ax4.plot([iteration[0], iteration[-1]], [pbounds[8][1], pbounds[8][1]], 'b--')
            ax5.plot([iteration[0], iteration[-1]], [pbounds[9][0], pbounds[9][0]], 'r--')
            ax5.plot([iteration[0], iteration[-1]], [pbounds[9][1], pbounds[9][1]], 'r--')
            ax5.plot([iteration[0], iteration[-1]], [pbounds[10][0], pbounds[10][0]], 'b--')
            ax5.plot([iteration[0], iteration[-1]], [pbounds[10][1], pbounds[10][1]], 'b--')
           

        # Put unit labels on axes where necessary. 
        ax1.set_ylabel(r"$R_E$", fontsize=8)
        ax2.set_ylabel(r"$ x10^{-5} eV cm^{-3} s^{-1}$", fontsize=8)
        
        if save: 
            
            fig.savefig(self.plot_path+"{}/{}_{}_model_parameter_changes_{}{}.png".format(self.current_model,self.filename, self.current_model, self.model['cost func'], savetag))

        self.fig_param = fig 

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
        xp_y, yp_y, zp_y, etad_y, xp_z, yp_z, zp_z, etad_z = gm.calculate_meridian_planes(self.model['x'], self.model['y'], self.model['z'], self.model['etad'])
        
        #Get meridian data for etam. 
        xp_y, yp_y, zp_y, etam_y, xp_z, yp_z, zp_z, etam_z = gm.calculate_meridian_planes(self.model['x'], self.model['y'], self.model['z'], self.model['etam'])
        
        # For the slice with constant y. 
        #y_uniq = abs(self.model['y'][0,:,0])
        #i_y = np.where(y_uniq == min(y_uniq))[0][0]
        # i_y = self.model.n[1]//2
        #xp_y = self.model['x'][:,i_y]
        #yp_y = self.model['y'][:,i_y]
        #zp_y = self.model['z'][:,i_y]
        #etad_y = self.model['etad'][:,i_y]
        #etam_y = self.model['etam'][:,i_y]
        #plane_value_y = self.model['y'][0,i_y,0]

        # For the slice with constant z. 
        #z_uniq = abs(self.model['z'][:,0,0])
        #i_z = np.where(z_uniq == min(z_uniq))[0][0]
        # i_z = self.model.n[2]//2
        #xp_z = self.model['x'][i_z]
        #yp_z = self.model['y'][i_z]
        #zp_z = self.model['z'][i_z]
        #etad_z = self.model['etad'][i_z]
        #etam_z = self.model['etam'][i_z]
        #plane_value_z = self.model['z'][i_z,0,0]
        
        # Calculate log10 eta values. If eta = 0, set log(eta) = vmin  
        letad_y = np.zeros(etad_y.shape)+vmin
        i = np.where(etad_y != 0)
        letad_y[i] = np.log10(etad_y[i])

        letam_y = np.zeros(etam_y.shape)+vmin
        i = np.where(etam_y != 0)
        letam_y[i] = np.log10(etam_y[i])

        letad_z = np.zeros(etad_z.shape)+vmin
        i = np.where(etad_z != 0)
        letad_z[i] = np.log10(etad_z[i])

        letam_z = np.zeros(etam_z.shape)+vmin
        i = np.where(etam_z != 0)
        letam_z[i] = np.log10(etam_z[i])
        
        
        # Now you can make the contour plot. 
        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(bottom=0.20, hspace=0.4, wspace=0.2)

        # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)

        # etad_y
        ax1 = fig.add_subplot(221)
        cont1 = ax1.contourf(xp_y, zp_y, letad_y, cmap=cmap, levels=levels, vmin=vmin, vmax=vmax)
        ax1.set_xlabel('X [RE]')
        ax1.set_ylabel('Z [RE]')
        ax1.set_title("n = {:.2f} cm".format(self.model['density'])+r"$^{-3}$"+"\nXZ Plane")
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
        cont2 = ax2.contourf(xp_y, zp_y, letam_y, cmap=cmap, levels=cont1.levels, vmin=vmin, vmax=vmax)
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
        cont3 = ax3.contourf(xp_z, yp_z, letad_z, cmap=cmap, levels=levels, vmin=vmin, vmax=vmax)
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
        cont4 = ax4.contourf(xp_z, yp_z, letam_z, cmap=cmap, levels=cont3.levels, vmin=vmin, vmax=vmax)
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
        for p,pval in enumerate(self.model['params best nm']):
                pv = pval 
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pv,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"

        fig.text(0.5, 0.02, label, ha='center')

        if save: 
            
            fig.savefig(self.plot_path+"{}/{}_data_{}_model_planes_opt{}.png".format(self.current_model,self.filename, self.current_model, savetag))
           
        self.fig = fig 

    def plot_earth_sun_line(self, save=False, savetag=""):
        '''This will plot the emissivity along a line close to the sun-earth line. 
        This will make it easier to see how the function compares to the simulation. 
        '''

        #Get Earth_sun line data for emissivity data. 
        xp, yp, zp, etad = gm.calculate_sunearth_line(self.model['x'], self.model['y'], self.model['z'], self.model['etad'])
		
		#Get Earth_sun line data for emissivity model. 
        xp, yp, zp, etam = gm.calculate_sunearth_line(self.model['x'], self.model['y'], self.model['z'], self.model['etam'])
        
        # For the slice with constant y. 
        #y_uniq = abs(self.model['y'][0,:,0])
        #i_y = np.where(y_uniq == min(y_uniq))[0][0]

        # For the slice with constant z. 
        #z_uniq = abs(self.model['z'][:,0,0])
        #i_z = np.where(z_uniq == min(z_uniq))[0][0]

        # Get data along sun-earth line. 
        #xp = self.model['x'][i_z,i_y]
        #yp = self.model['y'][i_z,i_y]
        #zp = self.model['z'][i_z,i_y]
        #etad = self.model['etad'][i_z,i_y]
        #etam = self.model['etam'][i_z,i_y]
        #plane_value_y = self.model['y'][0,i_y,0]
        #plane_value_z = self.model['z'][i_z,0,0]

        # Separate the model line into three colours for the different model sections. 
        # Separate the model line into three colours for the different model sections. 
        if self.current_model == "jorg":
        	i_msphere = np.where(xp <= self.model['params best nm'][0])
        	i_msheath = np.where((xp > self.model['params best nm'][0]) & (xp <= self.model['params best nm'][1]))
        	i_bow = np.where(xp > self.model['params best nm'][1])
        elif self.current_model == "cmem":
        	i_msphere = np.where(xp <= self.model['params best nm'][0]*self.model['r0lin'])
        	i_msheath = np.where((xp > self.model['params best nm'][0]*self.model['r0lin']) & (xp <= self.model['params best nm'][1]))
        	i_bow = np.where(xp > self.model['params best nm'][1])

        

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
        ax.set_title("Simulation Data vs {} Model - Sun-Earth Line\nn = {:.2f} cm".format(self.image_tag, self.model['density'])+r"$^{-3}$"+"\nOptimised Parameters")
        ax.set_xlim(0,25)

        # Add a label to show the model parameters. 
        label = ""
        info = gnau.get_parameter_info(model=self.current_model)
        parameter_names = [info[i][0] for i in info.keys()]
        parameter_units = [info[i][1] for i in info.keys()]
        for p,pval in enumerate(self.model['params best nm']):
                pv = pval 
                label += "{}={} {}, ".format(parameter_names[p], self.sig_figs(pv,3), parameter_units[p])
                if len(parameter_names)//2 == p+1:
                    label += "\n"
        fig.text(0.5, 0.02, label, ha='center')

        if save: 
            
            fig.savefig(self.plot_path+"{}/{}_data_{}_model_sunearth_opt{}.png".format(self.current_model,self.filename, self.current_model, savetag))
           


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
  
 
