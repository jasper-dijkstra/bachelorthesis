# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:49:40 2020

@author: Jasper Dijkstra

Initialize plume detection algorithm.

Adjust parameters to desired settings:
    - Target boundries: spatial extent of plume detection (lat: -90, 90; lon: -180, 180)
    - Target resolutions: grid cell resolution (km) of plume detection raster (> 7 km)
    - Basepath (folder containing inputs and where outputs will be stored)
    - Model behaviour settings:
        - Radius of buffer around TROPOMI plumes
        - Minimum amount of standard deviations within identification frames
        - Size of moving window frame
        - Steps between each moving window frame
        - Ignore all TROPOMI observations above sea
        - Incorporate CAMS (atmospheric composition) data in algorithm
        - Incorporate GFED and EDGAR data in the algorithm (to label plume origin)
    - GFED_path (folder where GFED file is stored)
    - EDGAR_path (folder where EDGAR file is stored)
    - CAMS_path (folder where CAMS file(s) are stored)
    

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


IMPORTANT NOTES:
    - Original filenames for GFED and EDGAR data needs to be retained
    - CAMS data needs to be stored using the format: 'EAC4_m[month]_d[day]_y[year].nc'
    - All input TROPOMI data needs to be stored as csv files in the folder: _basepath_\00_daily_csv.
        Use 'process_raw_data_to_daily_csv.py' to convert monthly csv files to daily required for the algorithm


"""


import os

import utilities as ut
import main as init

# ========================================================
# USER DEFINED PARAMETERS
# ========================================================

# == Boundary Conditions and resolutions == 
# Set the Target Boundaries (degrees)
lon_min = 100 #120 # 100 # minimum longitude
lon_max = 160 #126  # 160 # maximum longitude
lat_min = -50 #-24 # -50 # minimum latitude
lat_max = 0 #-18 # 0 # maximum latitude

# Setting the approximate target resolution (> 7)
lonres = 10 # km
latres = 10 # km

# == Model Behaviour Settings == 
buffersize = 4      # buffersize (radius of buffer around TROPOMI plumes)
stdevs = 1.8        # st.devs (minimum amount of st.devs within identification frames)
windowsize = 120    # windowsize (size of moving window frame in grid cells)
stepsize = 30       # stepsize (steps between each moving window frame in grid cells)

landsea = True      # Apply land-sea mask
cams = True         # Incorporate CAMS
gfededgar = True    # Compare results to GFED and EDGAR

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

daily_data = init.PlumeDetection(lat_min, lat_max, lon_min, lon_max, lonres, latres, \
                                 basepath, GFED_path, EDGAR_path, CAMS_path, \
                                 params = [buffersize, stdevs, windowsize, stepsize], \
                                     behaviour_settings = [landsea, cams, gfededgar, False, False])

