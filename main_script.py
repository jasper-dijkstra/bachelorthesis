# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:28:24 2020

@author: Jasper Dijkstra

This algorithm detects Carbon Monoxide plumes in TROPOMI data

Edit the parameters section to tune the algorithm as desired.

The rest of the script is built up as follows:
    - INITIALIZATION:
        TROPOMI data is processed to np.arrays
        If desired, measurements above oceans can be filtered out 
        Data is stored in dictionaries (daily_data_dict)
    - MAIN MODEL:
        Using a moving window, highest values are separated from background
        TROPOMI plumes can be compared to data from the Global Fire Emissions Database
    - HANDLING OUTPUT:
        If desired, textfiles containing plume coordinates are generated
        If desired, figures are generated

"""

#--------------------
# IMPORTING FUNCTIONS
#--------------------
import os
import copy
from datetime import datetime
import numpy as np

# Local imports
import handling_input as inpt
import handling_GFED as GFED
import handling_EDGAR as EDGAR
import handling_output_textfiles as txt
import handling_output_figures as figs
import handling_output as output
import masking_functions as mask
import raster_tools as raster
import fetching_winddata as wind
import utilities as ut



# ========================================================
# USER DEFINED PARAMETERS
# ========================================================

# == Boundary Conditions and resolutions == 
# Set the Target Boundaries (degrees)
lon_min = 100 #129 #120 #100 #100 #
lon_max = 160 #133 #126 #110 #160 #
lat_min = -50 #-22 #-24 #-30 #-50 #
lat_max = 0 #-18 #-40 #0 #
#boundaries = [-22, -18, 129, 133]
#boundaries = [-50, 0, 100, 160]

# Setting the approximate target resolution
lonres = 14 # km
latres = 14 # km

# == Desired Behaviour ==
# Apply operations:
use_wind_rotations = False   # rotate plumes in wind direction to improve results
apply_land_sea_mask = True  # filter out all TROPOMI data above the ocean
apply_overlap_filter = False # filter out plumes that occur for more than 3 consecutive satellite overpasses
    
# Outputs to be generated:
gen_txt_plume_coord = True  # txt file with plume coordinates
gen_figures = True          # figure with CO concentration (ppb) and mask (plume) as sublpot
gen_fig_wind_vector = False  # wind vector field figure


# == Directories ==
# Main Directory
basepath = ut.DefineAndCreateDirectory(r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\THESIS_WORKINGDIR\\')

# Directory where GFED files are stored
GFED_path = os.path.join(basepath, '01_GFED_hdf5_files' + os.path.sep) # Path to gfed .hdf5 files
EDGAR_path = os.path.join(basepath, '02_EDGAR_files' + os.path.sep) # Path to EDGAR .nc files

# ========================================================
# NON-USER DEFINED PARAMETERS
# ========================================================
# Store all defined extents in a list:
boundaries = [lat_min, lat_max, lon_min, lon_max]

# Calculating resolution based on lonres, latres
target_lon = int(abs((lon_max-lon_min)/(lonres/110)))
target_lat = int(abs((lat_max-lat_min)/(latres/110)))

# Output directories, for coordinates and figures
coord_dir = ut.DefineAndCreateDirectory(os.path.join(basepath + r'\04_output\plume_coordinates'))
fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath + r'\04_output\plume_figures'))

# Create a list with all files to apply the analysis on
input_files_directory = os.path.join(basepath + r'00_daily_csv\\')
files = ut.ListFilesInDirectory(input_files_directory, maxfiles=4)
del files[0:3] # Make sure script only runs for October 13, 2018 #DELETE THIS LATER

#%%
# ========================================================
# Initialization
# ========================================================
# Setting the time of starting the script
start = datetime.now()

# Reading daily csv files for specified area and day as np.arrays
daily_data = {}
for i, file in enumerate(files):    
    # Reading daily csv's as input array
    daily_data[i] = inpt.CSVtoArray(file, boundaries, target_lon, target_lat)
    
    # Filter measurements taken above the oceans (higher uncertainty)
    if apply_land_sea_mask == True:
        daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
        daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)

    # collect meteodata via ECMWF CDS API:
    if use_wind_rotations == True:
        u_wind, v_wind = wind.FetchWindData(daily_data[i], pressure=700, timerange=6, basepath=basepath)
        daily_data[i]['u_wind'] = u_wind
        daily_data[i]['v_wind'] = v_wind

print('Total time elapsed reading data: {}'.format(datetime.now()-start))



#%%
# ========================================================
# DETECTION ALGORITHM
# ========================================================
start_main = datetime.now()


""" Check for plumes based on three criteria: """
for day in daily_data:
    """ 1: At least 2 standard deviations above average of moving window """
    # Apply moving Window operation, with copy of CO_ppb
    arr = copy.deepcopy(daily_data[day]['CO_ppb'])
    plumes, co_average = raster.MovingWindow(arr, mask.identify_enhancements, window = (100,100), step = 20) 
    
    # Check if plume was detected in at least 95% of windows
    plumes[plumes >= 0.95] = 1 # If true, plume (1)
    plumes[plumes < 0.95] = 0 # If false, no plume (0)
    
    # Removing nuisances by making sure there is at least one neighbour
    neighbors = raster.CountNeighbors(plumes)  # Identify neighbors of each grid cell
    plumes[neighbors <= 1] = 0 # If there are 1 or fewer neighbors, undo identification as plume
    
    # Append to daily data dictionary
    #daily_data[day]['plume_mask'] = plumes
    
    """ 2: Check if plumes overlap > 3 days """
    if apply_overlap_filter:
        if day >= 2:
            # Identify plumes that occur > 3 days
            overlap = daily_data[day]['plume_mask'] * daily_data[day-1]['plume_mask'] \
                * daily_data[day-2]['plume_mask']
            valid = np.invert((overlap == 1)).astype(int)
            
            # Remove the values
            validity = lambda x: x * valid
            validity(daily_data[day-2]['plume_mask'])
            validity(daily_data[day-1]['plume_mask'])
            validity(daily_data[day]['plume_mask'])
               
        else:
            pass
    
    # PROBLEM: What if plume occurs 4 consecutive days? 0-3 are deleted,
    # next iteration 3 consecutive days are not detected
    

#for day in daily_data:
    """ 3: Check if plumes overlap with modelled GFED and EDGAR data, by drawing a buffer around TROPOMI data """
    
    # Read GFED and EDGAR data
    gfed_array = GFED.OpenGFED(GFED_path, boundaries, daily_data[day]['day'], \
                               daily_data[day]['month'], daily_data[day]['year'], lonres, latres)
    edgar_file = 'v432_CO_2012.0.1x0.1.nc' # EDGAR indsutry filename 
    edgar_array = EDGAR.OpenEDGAR(os.path.join(EDGAR_path + edgar_file), boundaries, lonres, latres)
    
    gfed_array[gfed_array > 0] = 10     # GFED plumes = 10
    edgar_array[edgar_array > 0] = 100  # EDGAR plumes = 100
    
    # Draw Buffer around TROPOMI plumes
    plumes_buffered = raster.DrawCircularBuffer(plumes, radius = 7) # Draw circular buffers around TROPOMI plumes

    # Delete GFED and EDGAR values outside of plume buffer
    gfed_tropomi = (plumes_buffered * gfed_array) # GFED within TROPOMI
    gfed_tropomi[gfed_tropomi > 0] = 1
    edgar_tropomi = (plumes_buffered * edgar_array) # EDGAR within TROPOMI
    edgar_tropomi[edgar_tropomi > 0] = 1
    
    # Adding all arrays with different plume origins:
    plumes = (plumes + gfed_array + edgar_array + gfed_tropomi + edgar_tropomi).astype(int)
    

    
    # Now: (0: no plume, 1: TROPOMI plume, 10: GFED plume (within bufferzone of TROPOMI plume, \
       # 11: TROPOMI + GFED identified plume))
    daily_data[day]['plume_mask'] = plumes

    
    """ 4: Rotate in the wind, to check if plume is detected downwind and not upwind """
    if use_wind_rotations:
        print('Wind rotations not available!')
        # Rotate identified plumes so they will be in line with the wind direction
        # If plume occurs exclusively downwind, it can be validated as plume.


print('Total time elapsed executing plume masking algorithm: {}'.format(datetime.now()-start_main))

#%%
# ========================================================
# GENERATING OUTPUTS
# ========================================================
start_end = datetime.now()

# Generate desired outputs based on earlier set preferences
for day in daily_data:
    if gen_txt_plume_coord:
        daily_data_dict = copy.deepcopy(daily_data[day])
        txt.NotePlumeCoordinates(daily_data_dict, coord_dir, lonres, latres)
    if gen_figures:
        figs.PlotFigures(daily_data[day], basepath, subplots=True)
    if gen_fig_wind_vector:
        assert use_wind_rotations == True, 'use_wind_rotations has to be True before wind vectorfield can be created!'
        output.CreateWindVector(daily_data[day], figtype = 'CO_ppb', figure_directory = fig_dir,\
                                labeltag = 'ppb', title='CO concentration and wind at 1000 hPa',\
                                    masking=True, skip=5)


print('Total time elapsed generating output: {}'.format(datetime.now()-start_end))
print('total time elapsed: {}'.format(datetime.now()-start))


