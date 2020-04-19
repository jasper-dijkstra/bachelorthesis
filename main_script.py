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
import os, sys
from datetime import datetime
import numpy as np

import read_write_functions as rw
import masking_functions as mask
import moving_window as window
import handling_output as output
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

# Decide whether or not a land-sea mask will be applied
gen_txt_plume_coord = False
gen_fig_xCO = False
gen_fig_plume = True

apply_land_sea_mask = True
apply_GEFD_data = True

# Setting the data working directory
basepath = r'C:\Users\jaspd\Desktop\THESIS_WORKINGDIR\\'

# Create a list with all files to apply the analysis on
input_csvs_dir = os.path.join(basepath + r'00_daily_csv\\')
files = []
for file in os.listdir(input_csvs_dir):
    if len(files) == 4: break # Limit the amount of input days
    if file.endswith(".csv"): files.append(os.path.join(input_csvs_dir + file))
    else:
        continue

#%%

#--------------------
# Initialization
#--------------------
# Setting the time of starting the script
start = datetime.now()

# Reading daily csv files for specified area and day as np.arrays
daily_data = {}
for i, file in enumerate(files):
        array = rw.reading_csv_as_nparray(file, boundaries, target_lon, target_lat)
        upd = {i : array}
        daily_data.update(upd)
        if apply_land_sea_mask == True:
            daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
            daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)
print('Finished reading data as np.array')

#%%
#--------------------
# Main Model
#--------------------


for day in daily_data:
    # Create a mask layer (1-0) for all enhancements, based on q-th percentile
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
