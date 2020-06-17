# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:49:40 2020

@author: Jasper Dijkstra

Initialize plume detection algorithm.

Adjust parameters to desired settings:
    - Target boundries (spatial extent of plume detection)
    - Target resolutions (grid cell resolution (km) of plume detection raster)
    - Basepath (folder where outputs i.a. will be stored)
    - GFED_path (folder where GFED file is stored (i.e. original filename should be retained))
    - EDGAR_path (folder where EDGAR file is stored (i.e. original filename should be retained))
    - !NOTE! all data to be used, need to be stored as csv files sorted per day in folder: 
        _basepath_\00_daily_csv


Outputs generated:
    - daily_data:
        python dictionary with all used data sorted per day
    - ASCII textfiles:
        text files of ASCII format, with all plume coordinates noted down
        saved at: _basepath_/04_output/plume_coordinates
    - Figures:
        Figures stored as Portable network graphics (png).
        saved at: _basepath_/04_output/plume_figures
        By default, one file is generated for each TROPOMI overpass.
        Each file contains two subplots:
            Total column CO as observed by TROPOMI (ppb)
            Enhanced data as identified by the algorithm (and GFED and EDGAR data for comparison)


Default behaviour of the model (NOTE this can be edited in the main.py script):
    - Incorporate wind data (False) [i.e. not working propoerly]
    - Apply overlap filter (filter plumes occuring consecutive days) (False) [i.e. not working propoerly]
    - Filter TROPOMI observations above oceans (True)
    - One output figure contains two subplots (CO concentration and plumes), these can also be separated in two different figures

"""


import os

import utilities as ut
import main as init

# ========================================================
# USER DEFINED PARAMETERS
# ========================================================

# == Boundary Conditions and resolutions == 
# Set the Target Boundaries (degrees)
lon_min = 100 #129 #120 #100 #100 # minimum longitude
lon_max = 160 #133 #126 #110 #160 # maximum longitude
lat_min = -50 #-22 #-24 #-30 #-50 # minimum latitude
lat_max = 0 #-18 #-40 #0 # maximum latitude

# Setting the approximate target resolution (> 7)
lonres = 10 # km
latres = 10 # km

# == Directories ==
# Main Directory
basepath = ut.DefineAndCreateDirectory(r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\THESIS_WORKINGDIR\\')

# Directory where GFED, EDGAR and CAMS files are stored
GFED_path = os.path.join(basepath, '01_GFED_hdf5_files' + os.path.sep) # Path to gfed .hdf5 files
EDGAR_path = os.path.join(basepath, '02_EDGAR_files' + os.path.sep) # Path to EDGAR .nc files
CAMS_path = os.path.join(basepath, '02_Store_CAMS_data' + os.path.sep) # Path to CAMS .nc files

# ========================================================
# START ALGORITHM
# ========================================================

daily_data = init.PlumeDetection(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath, GFED_path, EDGAR_path, CAMS_path)

