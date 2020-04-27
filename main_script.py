# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:28:24 2020

@author: Jasper Dijkstra

This script detects CO plumes in TROPOMI data

1. Reading data as 2D np.arrays
2. Separating enhancements in CO concentration from background
3. Rotate plumes (if possible) with wind direction, to find source
4. Plot figures/return coordinates of plume locations 


"""

#--------------------
# IMPORTING FUNCTIONS
#--------------------
import os
from datetime import datetime
import numpy as np

import handling_input as inpt
import handling_output as output
import masking_functions as mask
import moving_window as window
import utilities as ut


#--------------------
# PARAMETERS
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

# Decide what outputs have to be generated
gen_txt_plume_coord = False # txt file with plume coordinates
gen_fig_xCO = False # figure with CO concentration (ppb)
gen_fig_GFED = False # figure with GFED emissions (g C / m^2 / month)
gen_fig_plume = False # masked plume figure

gen_fig_GFED_buffer = False

# Decide whether or not land-sea mask and/or GFED data needs to be implemented
apply_land_sea_mask = True
compare_with_GFED = True

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
        # Should I add uncertainty of measurement?    
        day_data = inpt.reading_csv_as_nparray(file, boundaries, target_lon, target_lat)
        upd = {i : day_data}
        daily_data.update(upd)
        if apply_land_sea_mask == True:
            daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
            daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)

print('Total time elapsed reading data: {}'.format(datetime.now()-start))



#%%
#--------------------
# Main Model
#--------------------
start_main = datetime.now()

for day in daily_data:
    # Create a plume mask layer:
    arr = np.copy(daily_data[day]['CO_ppb'])
    outarr = window.MovingWindow(arr, window=(100,100), step=20, treshold=0.95)
    neighbors = window.CheckSurroundings(outarr) # Identify neighbrs of each grid cell
    outarr[neighbors < 1] = 0 # Removing nuisances by making sure there is at least one neighbour
    daily_data[day].update({'plume_mask':outarr, 'neighbors':neighbors})
    
    # Check if all plumes correspond with modelled GFED data
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
# Handling output
#--------------------
start_end = datetime.now()

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
    if gen_fig_GFED == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'GFED_emissions', title=None)
    if gen_fig_GFED_buffer == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'GFED_buffers', title=None)

print('Total time elapsed generating output: {}'.format(datetime.now()-start_end))
print('total time elapsed: {}'.format(datetime.now()-start))
