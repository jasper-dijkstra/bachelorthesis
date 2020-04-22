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
gen_fig_xCO = True # xCO figure
gen_fig_plume = True # masked plume figure

# Decide whether or not land-sea mask and/or GEFD data needs to be implemented
apply_land_sea_mask = False
apply_GFED_mask = True

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
        array = inpt.reading_csv_as_nparray(file, boundaries, target_lon, target_lat)
        upd = {i : array}
        daily_data.update(upd)
        if apply_land_sea_mask == True:
            daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
            daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)
        if apply_GFED_mask == True:
            daily_data[i]['CO_ppb'] = mask.GFED_mask(daily_data[i], 'CO_ppb')
            daily_data[i]['count_t'] = mask.GFED_mask(daily_data[i], 'count_t')

print('Total time elapsed reading data: {}'.format(datetime.now()-start))



#%%
#--------------------
# Main Model
#--------------------


for day in daily_data:
    # Create a plume mask layer:
    arr = np.copy(daily_data[day]['CO_ppb'])
    outarr = window.moving_window(arr, window=(100,100), step=20, treshold=0.95)
    daily_data[day].update({'plume_mask':outarr})
# Function(inputs = np.array per day plus looking x days in the past?)
        # Returns masklayer, 1 (enhanced) and 0 (background)

# Rectangle om plumemask heen, en daarbinnen roteren

# Now, use this in combination with wind rotation ->



#%%
#--------------------
# Handling output
#--------------------

for day in daily_data:
    if gen_txt_plume_coord == True:
        coord_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_coordinates'))
        output.NotePlumeCoordinates(daily_data[day], coord_dir)
    if gen_fig_xCO == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'xCO', title=None)
    if gen_fig_plume == True:
        fig_dir = ut.DefineAndCreateDirectory(os.path.join(basepath, r'plume_figures'))
        output.CreateFigue(daily_data[day], fig_dir, figtype = 'mask', title=None)

print('total time elapsed: {}'.format(datetime.now()-start))
