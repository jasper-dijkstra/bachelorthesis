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

# Local imports
import handling_input as inpt
import handling_output as output
import masking_functions as mask
import moving_window as window
import fetching_winddata as ERA5
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
use_wind_rotations = True   # rotate plumes in wind direction to improve results
apply_land_sea_mask = True  # filter out all TROPOMI data above the ocean
compare_with_GFED = True    # compare TROPOMI data with modelled wildfires from GFED
    
# Outputs to be generated:
gen_txt_plume_coord = False # txt file with plume coordinates
gen_fig_xCO = True          # figure with CO concentration (ppb)
gen_fig_plume = True        # masked plume figure

gen_fig_GFED = False         # figure with GFED emissions (g C / m^2 / month)
gen_fig_GFED_buffer = False # 

# Setting other parameters
max_unc_ppb = 50        # Maximum uncertainty in TROPOMI data

# Setting the data working directory
basepath = ut.DefineAndCreateDirectory(r'C:\Users\jaspd\Desktop\THESIS_WORKINGDIR')

# Create a list with all files to apply the analysis on
input_files_directory = os.path.join(basepath + r'00_daily_csv\\')
files = ut.ListCSVFilesInDirectory(input_files_directory, maxfiles=4)


#%%
#--------------------
# Initialization
#--------------------
# Setting the time of starting the script
start = datetime.now()

# Reading daily csv files for specified area and day as np.arrays
daily_data = {}
for i, file in enumerate(files):    
        day_data = inpt.reading_csv_as_nparray(file, boundaries, target_lon, target_lat, max_unc_ppb)
        upd = {i : day_data}
        daily_data.update(upd)
        if apply_land_sea_mask == True:
            daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
            daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)

# collect meteodata via ECMWF CDS API:
if use_wind_rotations == True:
    downloaded_nc_path = ERA5.DownloadERA5(daily_data, pressure_level=850)
    daily_data = ERA5.ProcessMeteo(daily_data, downloaded_nc_path)

print('Total time elapsed reading data: {}'.format(datetime.now()-start))



#%%
#--------------------
# MAIN MODEL
#--------------------
start_main = datetime.now()


# Check for plumes based on three criteria:
for day in daily_data:
    # 1: At least 2 standard deviations above average of moving window
    arr = np.copy(daily_data[day]['CO_ppb'])
    outarr = window.MovingWindow(arr, window=(100,100), step=20, treshold=0.95) # Apply moving Window operation
    neighbors = window.CheckSurroundings(outarr) # Identify neighbors of each grid cell
    outarr[neighbors < 1] = 0 # Removing nuisances by making sure there is at least one neighbour
    daily_data[day].update({'plume_mask':outarr, 'neighbors':neighbors})
    
    # 2: Rotate in the wind, to check if (center of) plumes overlap multiple days
        # Wildfires tend to occur less than one day!
    
    # 3: Check if all plumes overlap with a buffer around modelled GFED data
    if compare_with_GFED == True:
        gfed_array = inpt.ReadGFED(daily_data[day]) # Read GFED data as array
        gfed_buffered = window.DrawBuffer(gfed_array, buffersize = 10) # Draw circular buffers around GFED plumes
        outarr = daily_data[day]['plume_mask'] * gfed_buffered # Only allow plumes within buffersize
        daily_data[day].update({'plume_mask':outarr, 'GFED_emissions':gfed_array, 'GFED_buffers':gfed_buffered})

# Find center of plume (also, write this txt instead of every gridcell)
# Draw rectangle around it, size depending on size plume (find ideal size)
# Rotate within rectangle, so plume will be alligned with wind direction
# Reference this with plumes we found at other days???

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
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'CO_ppb', title=None)
    if gen_fig_plume == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'plume_mask', title=None)

# Is it really neccesary to plot these figures?
    if gen_fig_GFED == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'GFED_emissions', title=None)
    
    
    if gen_fig_GFED_buffer == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'GFED_buffers', title=None)

print('Total time elapsed generating output: {}'.format(datetime.now()-start_end))
print('total time elapsed: {}'.format(datetime.now()-start))
