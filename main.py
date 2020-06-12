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
import masking_functions as mask
import raster_tools as raster
import fetching_winddata as wind
import utilities as ut


# ========================================================
# NON-USER DEFINED PARAMETERS
# ========================================================

def InitializingParameters(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath):
    # Store all defined extents in a list:
    boundaries = [lat_min, lat_max, lon_min, lon_max]
    
    # Calculating resolution based on lonres, latres
    target_lon = int(abs((lon_max-lon_min)/(lonres/110)))
    target_lat = int(abs((lat_max-lat_min)/(latres/110)))
    
    # Create a list with all files to apply the analysis on
    input_files_directory = os.path.join(basepath + r'00_daily_csv\\')
    files = ut.ListFilesInDirectory(input_files_directory, maxfiles=4)
    del files[0:3] # Make sure script only runs for October 13, 2018 #DELETE THIS LATER

    return boundaries, target_lon, target_lat, files


# ========================================================
# Initialization (reading data)
# ========================================================

def CollectingData(boundaries, target_lon, target_lat, files, basepath, apply_land_sea_mask, use_wind_rotations):
    # Setting the time of starting the script
    start = datetime.now()
    
    # Reading daily csv files for specified area and day as np.arrays
    daily_data = {}
    for i, file in enumerate(files):    
        # Reading daily csv's as input array
        daily_data[i] = inpt.CSVtoArray(file, boundaries, target_lon, target_lat)
        
        # Filter measurements taken above the oceans (higher uncertainty)
        if apply_land_sea_mask:
            daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
            daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)
    
        # collect meteodata via ECMWF CDS API:
        if use_wind_rotations:
            u_wind, v_wind = wind.FetchWindData(daily_data[i], pressure=700, timerange=6, basepath=basepath)
            daily_data[i]['u_wind'] = u_wind
            daily_data[i]['v_wind'] = v_wind
    
    print('Total time elapsed reading data: {}'.format(datetime.now()-start))

    return daily_data


# ========================================================
# DETECTION ALGORITHM
# ========================================================

def Detection(params, daily_data, boundaries, GFED_path, EDGAR_path, lonres, latres, apply_overlap_filter, use_wind_rotations):
    start_main = datetime.now()
    
    
    """ Check for plumes based on three criteria: """
    for day in daily_data:
        """ 1: At least 2 standard deviations above average of moving window """
        # Apply moving Window operation, with copy of CO_ppb
        arr = copy.deepcopy(daily_data[day]['CO_ppb'])
        plumes, co_average = raster.MovingWindow(arr, mask.identify_enhancements, \
            window = (params[2],params[2]), step = params[3], st_devs = params[1]) 
        
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
        plumes_buffered = raster.DrawCircularBuffer(plumes, radius = params[0]) # Draw circular buffers around TROPOMI plumes
    
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
    
    return daily_data

# ========================================================
# GENERATING OUTPUTS
# ========================================================

def GeneratingOutputs(daily_data, basepath, lonres, latres, use_wind_rotations, gen_fig_wind_vector):
    start_end = datetime.now()
    
    # Generate desired outputs based on earlier set preferences
    for day in daily_data:
        # Generate text files with plume coordinates
        daily_data_dict = copy.deepcopy(daily_data[day])
        txt.NotePlumeCoordinates(daily_data_dict, basepath, lonres, latres)
        
        # Generate figures (subplots)
        figs.PlotFigures(daily_data[day], basepath, subplots=True)
        
        # Generate figures of the wind vectors
        if gen_fig_wind_vector:
            assert use_wind_rotations == True, 'use_wind_rotations has to be True before wind vectorfield can be created!'
            figs.CreateWindVector(daily_data[day], basepath, skip=30)
    
    
    print('Total time elapsed generating output: {}'.format(datetime.now()-start_end))
    
    return


def PlumeDetection(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath, GFED_path, EDGAR_path):
    start = datetime.now()
    print(f'Algorithm started at: {start}')
    print('')
    print('')
    
    # params = [
        # buffersize (radius of buffer around TROPOMI plumes),
        # st.devs (minimum amount of st.devs within identification frames),
        # windowsize (size of moving window frame in grid cells),
        # stepsize (steps between each moving window frame in grid cells)
        # ]
    params = [7, 1, 120, 20]
    
    apply_land_sea_mask = True
    apply_overlap_filter = False
    use_wind_rotations = False
    
    
    boundaries, target_lon, target_lat, files = InitializingParameters(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath)
    daily_data = CollectingData(boundaries, target_lon, target_lat, files, basepath, apply_land_sea_mask, use_wind_rotations)
    print('')
    print('')
    daily_data = Detection(params, daily_data, boundaries, GFED_path, EDGAR_path, lonres, latres, apply_overlap_filter, use_wind_rotations)
    print('')
    print('')
    GeneratingOutputs(daily_data, basepath, lonres, latres, use_wind_rotations, gen_fig_wind_vector=False,)
    
    print('')
    print('')
    print(f'Algorithm finished at: {start}')
    print(f'Total time elapsed: {datetime.now()-start}')
    print('')
    
    return daily_data