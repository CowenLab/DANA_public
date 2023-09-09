# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 18:41:32 2020

Iterates through all data directories. 
Relies on data collected by LD_Save_and_simplify_data_for_collaboration.m
That data can be found here...
https://www.dropbox.com/s/2sn8wzsg7blox8l/OK.zip?dl=0

For users other than Stephen, it's probably a good idea to just make a copy 
of this function with your own name at the end (e.g., LD_Session_Iterator_Susie.py) 
so that your changes to the paths will not interfer with the master version.

@author: Stephen Cowen
"""
import os
import glob
import sys
#import importlib
import pickle
# DIRS
CODE_DIR = r'C:\Users\Stephen Cowen\Documents\GitHub\Neural_Ensemble_Detection/Python/' # Yes, this is not the preferred python way to do this but it was the easiest until I understand Python a bit better.
PROCESSED_DATA_DIR = 'C:/Temp/OK' # download from https://www.dropbox.com/s/2sn8wzsg7blox8l/OK.zip?dl=0
ANA_RESULTS_DIR = 'C:/Temp/AnaResults' # Where you want all of the results from each experimental session/animal aggregated. Stored in .pkl files.

# Make a results directory if it does not exist.
if os.path.isdir(ANA_RESULTS_DIR) == False:
    os.mkdir(ANA_RESULTS_DIR)

# Add the target functions that we create to the python path. This will also allow me to import other packages in this directory.
sys.path.insert(0, CODE_DIR)
#from AutoEncoder_Demo import Analysis, Meta_Analysis
from Spike_by_Time_Matrices import Analysis, Meta_Analysis
#from Plot_Placefields import Analysis, Meta_Analysis
#from AutoEncode_Spike_Time import Analysis, Meta_Analysis
#importlib.reload(sys.modules['AutoEncoder_Demo'])
# Compile a list of all data directories.
animal_dirs = glob.glob(PROCESSED_DATA_DIR + "/8*")
all_dirs = []
for animal_dir in animal_dirs:
    ses_dirs = glob.glob(animal_dir + "/*")
    for ses_dir in ses_dirs:    
        all_dirs.append(ses_dir)

# Run the function in each directory (assumes the function cd's to that directory)
OUT = []
for a_dir in all_dirs:
    print (a_dir)
    tmp, day = os.path.split(a_dir)
    tmp, animal = os.path.split(tmp)
    fname = ANA_RESULTS_DIR + '/' + animal + '_' + day + '.pkl'
    # Analyze the data.
    OUT = Analysis(a_dir)
    # Save results to the output directory for meta-level analysis.
    f = open(fname,"wb")
    pickle.dump(OUT,f,pickle.HIGHEST_PROTOCOL)
    f.close()

# Run an analyses over all sessions and animals based on those pickle files stored above.
Meta_Analysis(ANA_RESULTS_DIR)