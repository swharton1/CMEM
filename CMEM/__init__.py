# This is the __init__.py file for the CMEM module. 

import os

# This is a backup to set the environment variables if they haven't been set externally. 
#This is mostly so I can do testing in ipython3. 
if "PLOT_PATH" not in os.environ: 
    os.environ["PLOT_PATH"] = "/home/s/sw682/Code/plots/CMEM_plots/"
if "PICKLE_PATH" not in os.environ:
    os.environ["PICKLE_PATH"] = "/data/sol-ionosphere/sw682/CMEM_pickled_models/" 
if "PPMLR_PATH" not in os.environ:
    os.environ["PPMLR_PATH"] = "/data/sol-ionosphere/SMILE/PPMLR/"
#os.environ["OPENGGCM_PATH"] = "/Users/sw682/Documents/Local_Code/CMEM/OpenGGCM/" 

from . import read_ppmlr
from . import visualise_nonopt
from . import fit_emissivity_models_old
from . import fit_emissivity_models 
from . import visualise_models
from . import compare_optimised_models 
from . import get_names_and_units 
from . import boundary_emissivity_functions 
from . import set_initial_params 
from . import get_meridians
from . import ppmlr_fits 
from . import coord_conv

 
