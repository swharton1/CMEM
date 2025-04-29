import numpy as np
import matplotlib.pyplot as plt
import os
import pickle


from SXI_Core import get_names_and_units as gnau 
from SXI_Core import get_meridians as gm 
from . import boundary_emissivity_functions as bef


class compare_models():
    def __init__(self):
        self.pickle_path = os.environ.get("PICKLE_PATH")
        self.plot_path = os.environ.get("PLOT_PATH")

    def read_pickle(self, model_tag, filename):
        '''This will read a single pickle file. '''

        with open(os.path.join(self.pickle_path, model_tag+"_optimised", filename), 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def compare_cost_between_models(self, save=True, savetag="", init_method=1, add_bad_fits=False, normalisation='normalised'):

        # List the filenames of the optimised models you want to compare. 
        if init_method == 2:
            filenames_jorg = ["S05D05V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D7.5V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D12.25V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D20V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D25V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D35V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",] 
        
            filenames_cmem = ["S05D05V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D7.5V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D12.25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D20V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D35V400B0000-05rad.fits_cmem_normalised_im2_.pkl",] 
        else:

            filenames_jorg = ["S05D05V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D7.5V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D12.25V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D20V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D25V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D35V400B0000-05rad.dat_jorg_normalised.pkl",] 
        
            filenames_cmem = ["S05D05V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D7.5V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D12.25V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D20V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D25V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D35V400B0000-05rad.fits_cmem_normalised_im1_.pkl",] 
        
        # Read these files. 
        pickle_dict_jorg = []
        for f in filenames_jorg:
            with open(os.path.join(self.pickle_path,"jorg_optimised",f), 'rb') as file:
                pickle_dict_jorg.append(pickle.load(file))

        pickle_dict_cmem = []
        for f in filenames_cmem:
            with open(os.path.join(self.pickle_path,"cmem_optimised",f), 'rb') as file:
                pickle_dict_cmem.append(pickle.load(file))
                       
        # Get minimum cost values out and density values out.
        density_jorg = np.array([p["density"] for p in pickle_dict_jorg])
        density_cmem = np.array([p["density"] for p in pickle_dict_cmem])

        cost_jorg = np.array([p["min cost"] for p in pickle_dict_jorg])
        cost_cmem = np.array([p["min cost"] for p in pickle_dict_cmem])

        # Here, we have added an extra step after comments from referee 2 about normalisation. 
        sum_eta_jorg = np.array([p["etad"].sum() for p in pickle_dict_jorg])
        sum_eta_cmem = np.array([p["etad"].sum() for p in pickle_dict_cmem])
        
        sum_squared_eta_jorg = np.array([(p["etad"]**2).sum() for p in pickle_dict_jorg])
        sum_squared_eta_cmem = np.array([(p["etad"]**2).sum() for p in pickle_dict_cmem])
        
        N_jorg = np.array([p["etad"].size for p in pickle_dict_jorg]) 
        N_cmem = np.array([p["etad"].size for p in pickle_dict_cmem]) 
        
        
        print (sum_eta_jorg)
        print (sum_eta_cmem) 
        print (sum_squared_eta_jorg)
        print (sum_squared_eta_cmem) 
        print (N_jorg)
        print (N_cmem) 
        
        #Apply chosen normalisation. 
        if normalisation == 'normalised':
            #This is the default I have been using to make the original plots. 
            cost_unit = ''
        elif normalisation == 'none': 
            #This removes all forms of normalisation. 
            cost_jorg = cost_jorg*sum_squared_eta_jorg
            cost_cmem = cost_cmem*sum_squared_eta_cmem 
            cost_unit = '[(eV cm'+r'$^{-3}$ s'+r'$^{-1}$)'+r'$^2$]'
        elif normalisation == 'eta':
            #This normalises by the sum of eta, instead of sum of eta squared. 
            cost_jorg = (cost_jorg*sum_squared_eta_jorg)/sum_eta_jorg
            cost_cmem = (cost_cmem*sum_squared_eta_cmem)/sum_eta_cmem 
            cost_unit = '[eV cm'+r'$^{-3}$ s'+r'$^{-1}$]'
        elif normalisation == 'N':
            #This just normalises by the number of data points, 
            #only taking into account the size of the cube. 
            cost_jorg = (cost_jorg*sum_squared_eta_jorg)/N_jorg
            cost_cmem = (cost_cmem*sum_squared_eta_cmem)/N_cmem 
            cost_unit = '[(eV cm'+r'$^{-3}$ s'+r'$^{-1}$)'+r'$^2$]'
        else:
            raise ValueError("Not picked a valid normalising constant: 'normalised', 'none' or 'eta'")    
            
            



        # Sort title and linestyle out. 
        if init_method == 2:
            title = "Variation of Cost with Solar Wind Density: Method 2"
            ls = "dashed"
           
        else:
            title = "Variation of Cost with Solar Wind Density: Method 1"
            ls = "dotted"
            

        # # Now create a plot. 
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(density_jorg, cost_jorg, "b", marker='x', linestyle=ls, label="Jorgensen")
        ax.plot(density_cmem, cost_cmem, "r", marker='x', linestyle=ls, label="CMEM")
        ax.set_xlabel("Density (cm"+r"$^{-3}$"+")")
        ax.set_ylabel("Minimum Cost "+cost_unit)
        ax.set_title(title)
        ax.legend(loc="best")
        ax.grid()

        # Add bad fits as black crosses. 
        if add_bad_fits:
            ax.plot(density_jorg[0], cost_jorg[0], marker='x', c='k', zorder=3)
            ax.plot(density_cmem[0], cost_cmem[0], marker='x', c='k', zorder=3)
            #ax.plot(density_cmem[1], cost_cmem[1], marker='x', c='k', zorder=3)

        #Add label to show normalisation method. 
        fig.text(0.9,0.03,normalisation, ha='right', fontsize=8) 

        if save:
            fig.savefig(self.plot_path+"cost_comparison{}_{}.png".format(savetag, normalisation))



    def compare_cost_between_models_both_methods(self, save=True, savetag="", add_bad_fits=False):

        # List the filenames of the optimised models you want to compare. 
       
        filenames_jorg_A1A2 = ["S05D05V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D7.5V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D12.25V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D20V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D25V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",
                          "S05D35V400B0000-05rad.dat_jorg_normalised_A1A2.pkl",] 
        
        filenames_cmem_A1A2 = ["S05D05V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D7.5V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D12.25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D20V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D35V400B0000-05rad.fits_cmem_normalised_im2_.pkl",] 
        

        filenames_jorg = ["S05D05V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D7.5V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D12.25V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D20V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D25V400B0000-05rad.dat_jorg_normalised.pkl",
                          "S05D35V400B0000-05rad.dat_jorg_normalised.pkl",] 
        
        filenames_cmem = ["S05D05V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D7.5V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D12.25V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D20V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D25V400B0000-05rad.fits_cmem_normalised_im1_.pkl",
                          "S05D35V400B0000-05rad.fits_cmem_normalised_im1_.pkl",] 
        
        # Read these files. 
        pickle_dict_jorg_A1A2 = []
        for f in filenames_jorg_A1A2:
            with open(os.path.join(self.pickle_path,"jorg_optimised",f), 'rb') as file:
                pickle_dict_jorg_A1A2.append(pickle.load(file))

        pickle_dict_cmem_A1A2 = []
        for f in filenames_cmem_A1A2:
            with open(os.path.join(self.pickle_path,"cmem_optimised",f), 'rb') as file:
                pickle_dict_cmem_A1A2.append(pickle.load(file))

        pickle_dict_jorg = []
        for f in filenames_jorg:
            with open(os.path.join(self.pickle_path,"jorg_optimised",f), 'rb') as file:
                pickle_dict_jorg.append(pickle.load(file))

        pickle_dict_cmem = []
        for f in filenames_cmem:
            with open(os.path.join(self.pickle_path,"cmem_optimised",f), 'rb') as file:
                pickle_dict_cmem.append(pickle.load(file))
                       
        # Get minimum cost values out and density values out.
        density_jorg_A1A2 = np.array([p["density"] for p in pickle_dict_jorg_A1A2])
        density_cmem_A1A2 = np.array([p["density"] for p in pickle_dict_cmem_A1A2])
        density_jorg = np.array([p["density"] for p in pickle_dict_jorg])
        density_cmem = np.array([p["density"] for p in pickle_dict_cmem])

        cost_jorg_A1A2 = np.array([p["min cost"] for p in pickle_dict_jorg_A1A2])
        cost_cmem_A1A2 = np.array([p["min cost"] for p in pickle_dict_cmem_A1A2])
        cost_jorg = np.array([p["min cost"] for p in pickle_dict_jorg])
        cost_cmem = np.array([p["min cost"] for p in pickle_dict_cmem])

        # Sort title out. 
        title = "Variation of Cost with Solar Wind Density: Both Methods"

        # # Now create a plot. 
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # if add_bad_fits:
        #     jcolor = ['k','b','b','b','b','b']
        #     ccolor = ['k','k','r','r','r','r']
        # else: 
        jcolor='b'
        ccolor='r'

        ax.plot(density_jorg_A1A2, cost_jorg_A1A2, color=jcolor, marker='x', linestyle="dashed", label="Jorgensen: M2", zorder=1)
        ax.plot(density_cmem_A1A2, cost_cmem_A1A2, color=ccolor, marker='x', linestyle="dashed", label="CMEM: M2", zorder=1)
        ax.plot(density_jorg, cost_jorg, "b", marker='x', linestyle="dotted", label="Jorgensen: M1", zorder=2)
        ax.plot(density_cmem, cost_cmem, "r", marker='x', linestyle="dotted", label="CMEM: M1", zorder=2)
        ax.set_xlabel("Density (cm"+r"$^{-3}$"+")")
        ax.set_ylabel("Minimum Cost")
        ax.set_title(title)
        ax.legend(loc="best")
        ax.grid()

        # Add bad fits as black crosses. 
        if add_bad_fits:
            ax.plot(density_jorg[0], cost_jorg[0], marker='x', c='k', zorder=3)
            ax.plot(density_cmem[0], cost_cmem[0], marker='x', c='k', zorder=3)
            #ax.plot(density_cmem[1], cost_cmem[1], marker='x', c='k', zorder=3)

        if save:
            fig.savefig(self.plot_path+"cost_comparison_both{}.png".format(savetag))







class parameter_relationships():
    #This create plots of the variation of four key parameters with solar wind density, as determined from the manually fitted files. 
    def __init__(self, model="jorg"):
        self.pickle_path = os.environ.get("PICKLE_PATH")
        self.plot_path = os.environ.get("PLOT_PATH")

        self.current_model = model
        
        # Get model tags. 
        if self.current_model == "jorg":
            #self.model_tag = "jorg"
            self.image_tag = "Jorg"
        elif self.current_model == "cmem":
            #self.model_tag = "cmem"
            self.image_tag = "CMEM"

    def read_pickle(self, model_tag, filename):
        '''This will read a single pickle file. '''

        with open(os.path.join(self.pickle_path, model_tag+"_manual", filename), 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict

    def plot_A1_A2_relationships(self, save=False, savetag=""):
        '''This will make a plot of A1 and A2 as a function of the solar wind density.
    
        Suggested way to run
        --------------------
        ppmlr_list = [ppmlr1, ppmlr2, ppmlr3, ppmlr4, ppmlr5, ppmlr6]
        params_list = [params0, params0_2, params0_3, params0_4, params0_5, params0_6]
        plot_A1_A2_relationships(ppmlr_list, params_list)

        '''

        # List the filenames of the manually optimised models you want to compare. 
        filenames_manual = ["S05D05V400B0000-05rad.dat_{}_manual.pkl".format(self.current_model),
                            "S05D7.5V400B0000-05rad.dat_{}_manual.pkl".format(self.current_model),
                            "S05D12.25V400B0000-05rad.dat_{}_manual.pkl".format(self.current_model),
                            "S05D20V400B0000-05rad.dat_{}_manual.pkl".format(self.current_model),
                            "S05D25V400B0000-05rad.dat_{}_manual.pkl".format(self.current_model),
                            "S05D35V400B0000-05rad.dat_{}_manual.pkl".format(self.current_model),
        ]
                          
        
        # Read these files. 
       

        pickle_dict = []
        for f in filenames_manual:
            with open(os.path.join(self.pickle_path,self.current_model+"_manual",f), 'rb') as file:
                pickle_dict.append(pickle.load(file))

        

        densities = []
        mp = []
        bs = []
        A1 = []
        A2 = []
        for p in pickle_dict:
            densities.append(p['density'])
            # mp is p0 in the Lin model. 
            mp.append(p["params0"][0])
            bs.append(p["params0"][1])
            A1.append(p["params0"][2])
            A2.append(p["params0"][3])
        densities = np.array(densities)
        mp = np.array(mp)
        bs = np.array(bs)
        A1 = np.array(A1)
        A2 = np.array(A2)


        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(left=0.2)
        ax1 = fig.add_subplot(211)
        if self.current_model == "cmem":
            ax1b = ax1.twinx()
            ax1b.plot(densities, mp, "b", marker="x", linestyle="dotted", label=r"$p_0$")
            ax1.plot(densities, bs, "r", marker="x", linestyle="dotted", label=r"${r_0}^{bs}$")
        else:
            ax1.plot(densities, mp, "b", marker="x", linestyle="dotted", label=r"${r_0}^{mp}$")
            ax1.plot(densities, bs, "r", marker="x", linestyle="dotted", label=r"${r_0}^{bs}$")
       
        ax2 = fig.add_subplot(212)
        ax2.plot(densities, A1, "b", marker="x", linestyle="dotted", label=r"$A_1$")
        ax2.plot(densities, A2, "r", marker="x", linestyle="dotted", label=r"$A_2$")
        
        ax1.set_ylabel(r"$R_E$")
        ax2.set_xlabel("Density (cm"+r"$^{-3}$"+")")
        ax2.set_ylabel(r"eV cm$^{-3}$ s$^{-1}$")
        ax1.set_title("Estimated Parameter Variation with Solar Wind Density\n{}".format(self.image_tag))
        ax1.grid()
        ax2.grid()

        # Add straight lines of best fit. 
        pmp = np.polyfit(densities, mp, 1)
        pbs = np.polyfit(densities, bs, 1)
        pA1 = np.polyfit(densities, A1, 1)
        pA2 = np.polyfit(densities, A2, 1)

        if self.current_model == "cmem":
            ax1b.plot(densities, pmp[0]*densities + pmp[1], "b", label="m = {:.4f}".format(pmp[0])+" c = {:.4f}".format(pmp[1]))
            ax1.plot(densities, pbs[0]*densities + pbs[1], "r", label="m = {:.2f}".format(pbs[0])+" c = {:.2f}".format(pbs[1]))
            ax1.legend(loc="right")
            ax1b.legend(loc="lower center")
        else:
            ax1.plot(densities, pmp[0]*densities + pmp[1], "b", label="m = {:.2f}".format(pmp[0])+" c = {:.2f}".format(pmp[1]))
            ax1.plot(densities, pbs[0]*densities + pbs[1], "r", label="m = {:.2f}".format(pbs[0])+" c = {:.2f}".format(pbs[1]))
            ax1.legend(loc="best")

        ax2.plot(densities, pA1[0]*densities + pA1[1], "b", label="m = {:.2f}".format(pA1[0]*1e5)+r"$x10^{-5},$"+" c = {:.2f}".format(pA1[1]*1e5)+r"$x10^{-5}$")
        ax2.plot(densities, pA2[0]*densities + pA2[1], "r", label="m = {:.2f}".format(pA2[0]*1e5)+r"$x10^{-5},$"+" c = {:.2f}".format(pA2[1]*1e5)+r"$x10^{-5}$")
        ax2.legend(loc="best")

        if save:
            fig.savefig(self.plot_path+"manual_param_variation_with_density_{}.png".format(self.current_model))
    
class optimal_parameter_relationships():
    '''This will plot the final optimised parameters as functions of solar wind density for the CMEM model.'''
    def __init__(self, model="cmem"):
        self.pickle_path = os.environ.get("PICKLE_PATH")
        self.plot_path = os.environ.get("PLOT_PATH")

        self.current_model = model


    def read_pickle(self, model_tag, filename):
        '''This will read a single pickle file. '''

        with open(os.path.join(self.pickle_path, model_tag+"_optimised", filename), 'rb') as f: 
            pickle_dict = pickle.load(f)
        return pickle_dict
    
    def plot_optimal_parameters(self, save=False):
        '''This will create a plot of all 12 parameters as functions of the solar wind density. 
        It uses the files generated by method 2 for initial parameter estimation. 
        It is just set up for the CMEM model at the moment. '''

         # List the filenames of the optimised models you want to compare. 
        filenames = ["S05D05V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D7.5V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D12.25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D20V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D35V400B0000-05rad.fits_cmem_normalised_im2_.pkl",] 
  
        # Read these files. 
        pickle_dict = []
        for f in filenames:
            with open(os.path.join(self.pickle_path,"cmem_optimised",f), 'rb') as file:
                pickle_dict.append(pickle.load(file))

        # Read out the parameters for each solar wind density. OPTIMAL 
        densities = []
        p0 = []
        bs = []
        A1 = []
        A2 = []
        B = []
        alpha = []
        beta = []
        p1 = []
        p2 = []
        p3 = []
        ay_bs = []
        az_bs = [] 

        for p in pickle_dict:
            densities.append(p['density'])
            # mp is p0 in the Lin model. 
            p0.append(p["params best nm"][0])
            bs.append(p["params best nm"][1])
            A1.append(p["params best nm"][2])
            A2.append(p["params best nm"][3])
            B.append(p["params best nm"][4])
            alpha.append(p["params best nm"][5])
            beta.append(p["params best nm"][6])
            p1.append(p["params best nm"][7])
            p2.append(p["params best nm"][8])
            p3.append(p["params best nm"][9])
            ay_bs.append(p["params best nm"][10])
            az_bs.append(p["params best nm"][11])


        param_list = [p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs]
        
        # Read out the parameters for each solar wind density. INITIAL 
        p0i = []
        bsi = []
        A1i = []
        A2i = []
        Bi = []
        alphai = []
        betai = []
        p1i = []
        p2i = []
        p3i = []
        ay_bsi = []
        az_bsi = [] 

        for p in pickle_dict:
            
            p0i.append(p["params0"][0])
            bsi.append(p["params0"][1])
            A1i.append(p["params0"][2])
            A2i.append(p["params0"][3])
            Bi.append(p["params0"][4])
            alphai.append(p["params0"][5])
            betai.append(p["params0"][6])
            p1i.append(p["params0"][7])
            p2i.append(p["params0"][8])
            p3i.append(p["params0"][9])
            ay_bsi.append(p["params0"][10])
            az_bsi.append(p["params0"][11])


        param_listi = [p0i, bsi, A1i, A2i, Bi, alphai, betai, p1i, p2i, p3i, ay_bsi, az_bsi]
        
        
        info = gnau.get_parameter_info(model=self.current_model)
        parameter_names = [info[i][0] for i in info.keys()]
        parameter_units = [info[i][1] for i in info.keys()]
        #param_symbols = [r"$p_0$", r"${r_0}^{bs}$", r"$A_1$", r"$A_2$", "B", r"$\alpha$", r"$\beta$", r"$p_1$", r"$p_2$", r"$p_3$", r"${\alpha_y}^{bs}$", r"${\alpha_z}^{bs}$"]
        #param_units = ["", r"$R_E$", r"$eV cm^{-3} s^{-1}$", r"$eV cm^{-3} s^{-1}$", "", "", "", "", "", "", "",""]
        letters = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"]
        # Now make the figure. 
        fig = plt.figure(figsize=(8,10))
        fig.subplots_adjust(left=0.15, hspace=0, wspace=0.1, top=0.93, right=0.85)

        # Loop through to make the axes. 
        for p, pval in enumerate(param_list):
            ax = fig.add_subplot(6,2,p+1)
            if (p == 2) or (p == 3):
                ax.plot(densities, np.array(pval)*1e5, marker='x', linestyle="dotted", color="b")
                ax.plot(densities, np.array(param_listi[p])*1e5, marker='x', linestyle="dotted", color="r")
                ax.set_ylabel(parameter_names[p]+" "+r"$ (x10^{-5})$ "+parameter_units[p],fontsize=8)
            else:
                ax.plot(densities, pval, marker='x', linestyle="dotted", color="b")
                ax.plot(densities, param_listi[p], marker='x', linestyle="dotted", color="r")
                ax.set_ylabel(parameter_names[p]+" "+parameter_units[p])
            ax.set_xticks([0,5,10,15,20,25,30,35,40])
            ax.set_xlim(0,40)
            ax.grid()
            ax.text(0.02, 0.98, letters[p], transform=ax.transAxes, ha="left", va="top")
            if p >= 10: 
                ax.set_xlabel("Density (cm"+r"$^{-3}$"+")", fontsize=8)
            else:
                ax.xaxis.set_tick_params(labelbottom=False)
            
            if p%2 == 1: 
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")

            for label in (ax.get_xticklabels() + ax.get_yticklabels()): 
                label.set_fontsize(8)

            fig.text(0.5, 0.95, "Parameter Variation for the CMEM Model", ha="center")
            fig.text(0.85, 0.95, 'Initial', ha="right", color="r", fontsize=8)
            fig.text(0.85, 0.94, 'Optimised', ha="right", color="b", fontsize=8)

        if save:
            fig.savefig(self.plot_path+"optimal_param_variation_with_density_CMEM.png")
    
class magnetopause_model():
    '''This will create an Earth-Sun line plot, or extract the data from it, and work out
    the magnetopause position. Can use several definitions. '''

    def __init__(self):
        '''This takes in the threed model object you've created from a pickle file. '''
        
        self.pickle_path = os.environ.get("PICKLE_PATH")
        self.plot_path = os.environ.get("PLOT_PATH")

        #Get the names of the variables and units for plotting. 
        info = gnau.get_parameter_info(model='cmem')
        self.parameter_names = [info[i][0] for i in info.keys()]
        self.parameter_units = [info[i][1] for i in info.keys()]
        
        self.get_pickle_files()
        self.get_earth_sun_line_data()
        self.get_magnetopause_positions()
        
       

    def __repr__(self):
        return f"magnetopause model object."
    
    
    def get_pickle_files(self):
        '''This will read in the relevant pickle file for each simulation from the CMEM model. '''

        self.filenames = ["S05D05V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D7.5V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D12.25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D20V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D25V400B0000-05rad.fits_cmem_normalised_im2_.pkl",
                          "S05D35V400B0000-05rad.fits_cmem_normalised_im2_.pkl",] 

        # Read these files. 
        self.pickle_dict = []
        for f in self.filenames:
            with open(os.path.join(self.pickle_path,"cmem_optimised",f), 'rb') as file:
                self.pickle_dict.append(pickle.load(file))

    def get_earth_sun_line_data(self):
        '''This will calculate the important emissivity data along the Earth-Sun line. '''

        self.xp_array = []
        self.etam_array = []
        self.etad_array = []
        self.bs = []
        self.params = []
        self.density = []
        self.pdyn = []
        self.pmag = []
        self.bz = [] 

        for p, pval in enumerate(self.pickle_dict):
            
            #Get Earth_sun line data for emissivity data. 
            xp, yp, zp, etad = gm.calculate_sunearth_line(pval['x'], pval['y'], pval['z'], pval['etad'])
        
            #Get Earth_sun line data for emissivity model. 
            xp, yp, zp, etam = gm.calculate_sunearth_line(pval['x'], pval['y'], pval['z'], pval['etam'])

            # Get bowshock parameter too. You will need it. 
            bs = pval["params best nm"][1]
            params = pval["params best nm"]

            # Get density. 
            density = pval["density"]
            pdyn = pval["pdyn"]
            pmag = pval["pmag"]
            #bz = pval["bz"] 
            
            # Get parameter label information for plotting. 
            #if p == 0: 
            #    self.param_list = pval["parameter list"]
            #    self.param_unit = pval["parameter unit list"]

            self.xp_array.append(xp)
            self.etam_array.append(etam)
            self.etad_array.append(etad)
            self.bs.append(bs)
            self.params.append(params)
            self.density.append(density)
            self.pdyn.append(pdyn)
            self.pmag.append(pmag)
            #self.bz.append(bz) 

    def get_magnetopause_positions(self):
        '''This will go through the Earth-Sun line data and work out the magnetopause position from
        the data for all four techniques. '''

        # Set empty arrays. 
        self.r_cmem_array = []
        self.maxIx_array = []
        self.maxdIx_array = []
        self.f_array = []

        # Loop through each simulation. 
        for p, pval in enumerate(self.xp_array):

            xp = pval 
            etam = self.etam_array[p]
            etad = self.etad_array[p]

            # Get CMEM magnetopause first. OLD DEFINITION. 
            #xp_index = np.where(etam == etam.max())
            #self.r_cmem_array.append(xp[xp_index])
            #Get Lin coefficients. 
            
            lin_coeffs = bef.get_lin_coeffs(0, self.pdyn[p], self.pmag[p], -5.0)
            #Get Lin magnetopause for a range of theta and phi 
            rmp = bef.lin_scaled_func(0, 0, *lin_coeffs, p0=self.params[p][0], p1=self.params[p][7], p2=self.params[p][8], p3=self.params[p][9]) 
            self.r_cmem_array.append(rmp) 
            
            # Get max Ix second. 
            ix_index = np.where(etad == etad.max())
            self.maxIx_array.append(xp[ix_index])

            # Get max dIx third. 
            # Get difference between etad values. 
            dIx = np.array([etad[i+1] - etad[i] for i in range(len(etad) - 1)])

            # Get centre point positions for radial direction. 
            xp_cent = xp + (xp[1]-xp[0])/2
            xp_cent = xp_cent[0:-1]

            dix_index = np.where(dIx == dIx.max())
            self.maxdIx_array.append(xp[dix_index])

            # Get f=0.25 
            dr = xp[ix_index] - xp[dix_index]
            f = xp[dix_index] + 0.25*dr 
            self.f_array.append(f)

        self.r_cmem_array = np.array(self.r_cmem_array)
        self.maxIx_array = np.array(self.maxIx_array)
        self.maxdIx_array = np.array(self.maxdIx_array)
        self.f_array = np.array(self.f_array)
        
            



    def plot_earth_sun_line(self, sim=0, xlim=[0,25], save=False):
        '''This will read in the data for all six 
        
        Parameters
        ----------
        sim - simulation number index (0,1,2,3,4,5)
        xlim - def = [0,25]
        '''

        # Select data from the right simulation. 
        xp = self.xp_array[sim]
        etam = self.etam_array[sim]
        etad = self.etad_array[sim]

        
        # Select sections to colour. 
        i_msphere = np.where(xp < self.r_cmem_array[sim])
        i_msheath = np.where((xp >= self.r_cmem_array[sim]) & (xp < self.bs[sim]))
        i_bow = np.where(xp >= self.bs[sim])

         # Now you can make the line plot. 
        fig = plt.figure(figsize=(8,6))
        fig.subplots_adjust(bottom=0.20)
        ax = fig.add_subplot(111)

        ax.plot(xp, etad, 'k')
        ax.plot(xp[i_msphere], etam[i_msphere], 'b')
        ax.plot(xp[i_msheath], etam[i_msheath], 'r')
        ax.plot(xp[i_bow], etam[i_bow], 'g')

        ax.set_xlabel('X [RE]')
        ax.set_ylabel(r"eV cm$^{-3}$ s$^{-1}$")
        
        ax.set_title("Simulation Data vs CMEM Model - Sun-Earth Line\nn = {} cm".format(self.density[sim])+r"$^{-3}$"+"\n Optimised Parameters")
        
        ax.set_xlim(xlim)

        label = ""
        #info = gnau.get_parameter_info(model='cmem')
        #parameter_names = [info[i][0] for i in info.keys()]
        #parameter_units = [info[i][1] for i in info.keys()]
        for p,pval in enumerate(self.params[sim]):
                pv = pval 
                label += "{}={} {}, ".format(self.parameter_names[p], self.sig_figs(pv,3), self.parameter_units[p])
                if len(self.parameter_names)//2 == p+1:
                    label += "\n"
                    


        fig.text(0.5, 0.02, label, ha='center')

        # Now add on the magnetopause boundaries. 
        ax.plot([self.r_cmem_array[sim], self.r_cmem_array[sim]], [0, self.etam_array[sim].max()], 'r--', label=r"$r_{CMEM}$ = "+"{:.3f}".format(self.r_cmem_array[sim]))
        ax.text(self.r_cmem_array[sim], self.etam_array[sim].max(), r"$r_{CMEM}$", va="bottom", ha='right')

        ax.plot([self.maxIx_array[sim], self.maxIx_array[sim]], [0, self.etad_array[sim].max()], 'k--', label=r"max $Ix$ = "+"{:.3f}".format(self.maxIx_array[sim][0]))
        ax.text(self.maxIx_array[sim], self.etad_array[sim].max(), r"max $Ix$", va="bottom", ha='center')

        ax.plot([self.maxdIx_array[sim], self.maxdIx_array[sim]], [0, self.etad_array[sim].max()/2], 'b--', label=r"max d$Ix$ = "+"{:.3f}".format(self.maxdIx_array[sim][0]))
        ax.text(self.maxdIx_array[sim], self.etad_array[sim].max()/2, r"max d$Ix$", va="bottom", ha='right')

        ax.plot([self.f_array[sim], self.f_array[sim]], [0, self.etad_array[sim].max()*0.75], 'g--', label=r"$f_{0.25}$ = "+"{:.3f}".format(self.f_array[sim][0]))
        ax.text(self.f_array[sim], self.etad_array[sim].max()*0.75, r"$f_{0.25}$", va="bottom", ha='right')

        ax.legend()

        if save: 
            fig.savefig(self.plot_path+"earthsun_with_magnetopause_{}.png".format(sim))
        
    def plot_all_boundaries(self, save=False):
        '''This will plot r_cmem and the other magnetopause boundary definitions for each density. '''

        # Make the figure. 
        fig = plt.figure(figsize=(6,8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        ax1.plot(self.density, self.r_cmem_array, c='r', linestyle="dotted", marker='x', label=r"$r_{CMEM}$")
        ax1.plot(self.density, self.maxIx_array, c='k', linestyle="dotted", marker='^', label=r"max $Ix$")
        ax1.plot(self.density, self.f_array, c='g', linestyle="dotted", marker='s', label=r"$f_{0.25}$")
        ax1.plot(self.density, self.maxdIx_array, c='b', linestyle="dotted", marker='o', label=r"max d$Ix$")
        
        # ax1.set_xlabel("Density (cm"+r"$^{-3}$"+")")
        ax1.set_xlim(0,40)
        ax1.set_ylabel(r"$R_E$")
        ax1.legend(fontsize=8)
        ax1.set_title("Subsolar Magnetopause Boundary Positions", fontsize=12)
        ax1.grid()
        ax1.text(0.02, 0.97, "(a)", transform=ax1.transAxes, ha="left", va="top")
        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()): 
                label.set_fontsize(8)

        # Calculate mean differences. 
        self.mean_r = np.mean(self.r_cmem_array-self.r_cmem_array)
        self.mean_maxIx = np.mean(self.maxIx_array[:,0]-self.r_cmem_array)
        self.mean_maxdIx = np.mean(self.maxdIx_array[:,0]-self.r_cmem_array)
        self.mean_f = np.mean(self.f_array[:,0]-self.r_cmem_array)

        ax2.plot(self.density, self.r_cmem_array-self.r_cmem_array, c='r', linestyle="dotted", marker='x', label="Mean = {:.3f}".format(self.mean_r)+r"$R_E$")
        ax2.plot(self.density, self.maxIx_array[:,0]-self.r_cmem_array, c='k', linestyle="dotted", marker='^', label="Mean = {:.3f}".format(self.mean_maxIx)+r"$R_E$")
        ax2.plot(self.density, self.f_array[:,0]-self.r_cmem_array, c='g', linestyle="dotted", marker='s', label="Mean = {:.3f}".format(self.mean_f)+r"$R_E$")
        ax2.plot(self.density, self.maxdIx_array[:,0]-self.r_cmem_array, c='b', linestyle="dotted", marker='o', label="Mean = {:.3f}".format(self.mean_maxdIx)+r"$R_E$")
        ax2.set_xlabel("Density (cm"+r"$^{-3}$"+")")
        ax2.set_xlim(0,40)
        ax2.set_ylabel(r"$\Delta R_E$")
        ax2.set_ylim(-0.8, 1.0)
        ax2.legend(fontsize=8)
        ax2.set_title("Relative Subsolar Magnetopause Boundary Positions", fontsize=12)
        ax2.grid()
        ax2.text(0.02, 0.97, "(b)", transform=ax2.transAxes, ha="left", va="top")
        for label in (ax2.get_xticklabels() + ax2.get_yticklabels()): 
                label.set_fontsize(8)

        if save: 
            fig.savefig(self.plot_path+"magnetopause_positions.png")

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
