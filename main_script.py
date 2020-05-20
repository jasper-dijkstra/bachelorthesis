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
        !wind rotations!
        If desired, TROPOMI plumes can be compared to data from the Global Fire Emissions Database
    - HANDLING OUTPUT:
        If desired, textfiles containing plume coordinates are generated
        If desired, figures are generated

"""


#--------------------
# IMPORTING FUNCTIONS
#--------------------
import os
from datetime import datetime
import numpy as np
import scipy as sp

# Local imports
import handling_input as inpt
import handling_output as output
import masking_functions as mask
import raster_tools as raster
import fetching_winddata as wind
import utilities as ut

#--------------------
# PARAMETERS (user input)
#--------------------
# Set the Target Boundaries (degrees)
lon_min = 100
lon_max = 160
lat_min = -50
lat_max = 0
boundaries = [lat_min, lat_max, lon_min, lon_max]

# Setting the target resolution to ~7x7km
target_lon = int(abs((lon_max-lon_min)/(7/110)))
target_lat = int(abs((lat_max-lat_min)/(7/110)))

# Apply operations:
use_wind_rotations = False   # rotate plumes in wind direction to improve results
apply_land_sea_mask = True  # filter out all TROPOMI data above the ocean
    
# Outputs to be generated:
gen_txt_plume_coord = False  # txt file with plume coordinates
gen_fig_xCO = True          # figure with CO concentration (ppb)
gen_fig_plume = True        # masked plume figure
gen_fig_wind_vector = False  # wind vector field figure

# Setting other parameters
max_unc_ppb = 50        # Maximum uncertainty in TROPOMI data

# Setting the data working directory
basepath = ut.DefineAndCreateDirectory(r'C:\Users\jaspd\Desktop\THESIS_WORKINGDIR')

# Create a list with all files to apply the analysis on
input_files_directory = os.path.join(basepath + r'00_daily_csv\\')
files = ut.ListCSVFilesInDirectory(input_files_directory, maxfiles=4)
del files[0:3]

#%%
#--------------------
# Initialization
#--------------------
# Setting the time of starting the script
start = datetime.now()

# Reading daily csv files for specified area and day as np.arrays
daily_data = {}
for i, file in enumerate(files):    
    # Reading daily csv's as input array
    daily_data[i] = inpt.CSVtoArray(file, boundaries, target_lon, target_lat, max_unc_ppb)
    
    # Filter measurements taken above the oceans (higher uncertainty)
    if apply_land_sea_mask == True:
        daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
        daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)

    # collect meteodata via ECMWF CDS API:
    if use_wind_rotations == True:
        u_wind, v_wind = wind.FetchWindData(daily_data[i], pressure=1000, timerange=5)
        daily_data[i]['u_wind'] = u_wind
        daily_data[i]['v_wind'] = v_wind

print('Total time elapsed reading data: {}'.format(datetime.now()-start))



#%%
#--------------------
# MAIN MODEL
#--------------------
start_main = datetime.now()


""" Check for plumes based on three criteria: """
for day in daily_data:
    """ 1: At least 2 standard deviations above average of moving window """
    # Apply moving Window operation, with copy of CO_ppb
    arr = np.copy(daily_data[day]['CO_ppb'])
    plumes = raster.MovingWindow(arr, mask.identify_enhancements_3, window = (100,100), step = 20) 
    
    # Check if plume was detected in at least 95% of windows
    plumes[plumes >= 0.95] = 1 # If true, plume (1)
    plumes[plumes < 0.95] = 0 # If false, no plume (0)
    
    # Removing nuisances by making sure there is at least one neighbour
    neighbors = raster.CountNeighbors(plumes)  # Identify neighbors of each grid cell
    plumes[neighbors <= 1] = 0 # If there are 1 or fewer neighbors, undo identification as plume
    
    """ 2: Rotate in the wind, to check if (center of) plumes overlap multiple days """
    if use_wind_rotations:
        print('Wind rotations under construction!')
        # Wildfires tend to occur less than one day!
        # Draw rectangle around it, size depending on size plume (find ideal size)
        # Rotate within rectangle, so plume will be alligned with wind direction
        # Reference this with plumes we found at other days???

        
 
    """ 3: Check if plumes overlap with modelled GFED data, by drawing a buffer around TROPOMI data """

    gfed_array = inpt.ReadGFED(daily_data[day]) # Read GFED data as array
    
    # Buffer TROPOMI 
    plumes_buffered = raster.DrawCircularBuffer(plumes, radius = 7) # Draw circular buffers around TROPOMI plumes
    GFED_plumes = plumes_buffered * gfed_array # Array with overlapping GFED and TROPOMI plumes
    GFED_plumes[GFED_plumes > 0] = 10 # Setting all GFED plumes to value 10
    
    
# =============================================================================
#     # Buffer GFED 
#     gfed_buffered = raster.DrawCircularBuffer(gfed_array, radius = 7) # Draw circular buffers around TROPOMI plumes
#     GFED_plumes = gfed_buffered * plumes # Array with overlapping GFED and TROPOMI plumes
#     GFED_plumes[GFED_plumes > 0] = 10 # Setting all GFED plumes to value 10
# =============================================================================
    
    
    
    # Now do the same for steel data?
    
    # Adding all arrays with different plume origins:
    plumes = plumes + GFED_plumes
    
    # Now: (0: no plume, 1: TROPOMI plume, 10: GFED plume (within bufferzone of TROPOMI plume, \
       # 11: TROPOMI + GFED identified plume))
    
    
    daily_data[day]['plume_mask'] = plumes


print('Total time elapsed executing plume masking algorithm: {}'.format(datetime.now()-start_main))

#%%
#--------------------
# HANDLING OUTPUT
#--------------------
start_end = datetime.now()

# Generate desired outputs based on earlier set preferences
for day in daily_data:
    if gen_txt_plume_coord == True:
        coord_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_coordinates'))
        output.NotePlumeCoordinates(daily_data[day], coord_dir)
    if gen_fig_xCO == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateColorMap(daily_data[day], 'CO_ppb', fig_dir, labeltag = 'ppb')
    if gen_fig_plume == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateMaskMap(daily_data[day], 'plume_mask', fig_dir)
    if gen_fig_wind_vector == True:
        assert use_wind_rotations == True, 'use_wind_rotations has to be True before wind vectorfield can be created!'
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateWindVector(daily_data[day], figtype = 'CO_ppb', figure_directory = fig_dir,\
                                labeltag = 'ppb', title='CO concentration and wind at 1000 hPa',\
                                    masking=True, skip=30)



print('Total time elapsed generating output: {}'.format(datetime.now()-start_end))
print('total time elapsed: {}'.format(datetime.now()-start))


