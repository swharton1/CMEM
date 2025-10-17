# This is the __init__.py file for the CMEM module. 

import os

# This is a backup to set the environment variables if they haven't been set externally. 
#This is mostly so I can do testing in ipython3. 
if "PLOT_PATH" not in os.environ: 
    os.environ["PLOT_PATH"] = "/home/s/sw682/Code/plots/CMEM_plots/"
if "PICKLE_PATH" not in os.environ:
    os.environ["PICKLE_PATH"] = "/data/sol-ionosphere/sw682/pickled_files/CMEM_pickled_models/" 
if "PPMLR_PATH" not in os.environ:
    os.environ["PPMLR_PATH"] = "/data/smile/PPMLR/"


from . import visualise_nonopt
from . import fit_emissivity_models_old #Redundant 
from . import fit_emissivity_models 
from . import visualise_models
from . import compare_optimised_models 
from . import boundary_emissivity_functions 
from . import SXI_Core



 
