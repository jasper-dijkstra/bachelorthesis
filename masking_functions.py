# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:30:07 2020

@author: Jasper Dijkstra

This script contains functions that return as masks over TROPOMI CO data (or any 2d np.array)

Functions:
    1. land/sea mask
    2. Carbon mask (TO BE MADE)
    3. Plume mask (3 functions)
    
"""

import numpy as np
from global_land_mask import globe
from GFED_fire_emissions_mask import hdf5_to_mask as gfed

def land_sea_mask(array, boundaries):
    """
    
    Parameters
    ----------
    array : np.array
        geolocated 2d np.array that describes the values to be tested for land/sea.
    boundaries : list
        list defining the boundaries of array: [lat_min, lat_max, lon_min, lon_max].

    Returns
    -------
    None.

    """
    # Defining the lat and longitudinal ranges
    latrange = np.linspace(boundaries[0], boundaries[1], len(array))
    lonrange = np.linspace(boundaries[2], boundaries[3], len(array[0]))
    
    # Check for all lat/lon values if they are land (1) or water (0)
    mask = []
    for lat in range(len(latrange)):
        is_on_land = globe.is_land(latrange[lat], lonrange)
        is_on_land = is_on_land.astype(int)
        mask.append(is_on_land)
    mask = np.array(mask)

    land_arr = mask * array
    
    return land_arr


def GFED_mask(daily_data_dict, array_name):
    
    # Make sure array_name is correct
    supported_arrays = ['CO_ppb', 'count_t']
    if array_name not in supported_arrays:
        print('array_name was not recognised by GFED_mask()')
        # Write to logging file
        return
    
    # Define lat/lon ranges of target
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Define timescope
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    
    # Defining the lat and longitudinal ranges
    latrange = np.linspace(lat_min, lat_max, len(daily_data_dict['CO_ppb']))
    lonrange = np.linspace(lon_min, lon_max, len(daily_data_dict['CO_ppb'][0]))
    
    # Check for all lat/lon values if they correspond with wildfire emissions (1) or not (0)
    mask = []
    array = daily_data_dict[array_name]
    for lat in range(len(latrange)):
        #gfed.read_hdf5(day, month, year)
        is_wildfire = gfed.is_wildfire(latrange[lat], lonrange, day, month, year)
        is_wildfire = is_wildfire.astype(int)
        mask.append(is_wildfire)
    mask = np.array(mask)

    mask = np.multiply(mask, array)
    
    return mask



def identify_enhancements(frame_array, q):
    """
    
    Parameters
    ----------
    arr : np.array
        2D np.array.
    q : float (0-1)
        Quantile that has to be considered as peak.

    Returns
    -------
    arr : np.array
        2D np.array, where the background has been removed (lower quantile) and the enhancements have been identified (upper quantile)

    """
    frame_array[frame_array == 0] = 'nan' # Change all zeros to 'nan' so they won't be taken into account
    background = np.nanquantile(frame_array,q) # Define background concentration
    frame_array = frame_array - background # Remove background from actual concentation
    frame_array = np.nan_to_num(frame_array) # Rechange all 'nan' values to 0
    frame_array[frame_array < 0] = 0 # Change all negative values to 0
    frame_array[frame_array > 0] = 1 # Change all enhancements to 1
    return frame_array


def identify_enhancements_2(arr, q, st_devs = 1):
    arr[arr == 0] = 'nan' # Change all zeros to 'nan' so they won't be taken into account
    
    #average_whole_frame = np.nanmean(arr)
    stdev = st_devs*np.nanstd(arr)
    
    background = np.nanquantile(arr,q) # Define background concentration
    
    arr = arr - background # Remove background from actual concentation
    
    arr = arr - stdev # Remove also standard deviation, to make sure it really are peaks
    
    #arr[arr < 0] = 'nan' # Change all negative values to 0
    
    #average_quantile_frame = np.nanmean(arr)
    
    arr = np.nan_to_num(arr) # Rechange all 'nan' values to 0
    arr[arr < 0] = 0 # Change all negative values to 0
    arr[arr > 0] = 1 # Change all enhancements to 1
    return arr

def identify_enhancements_3(arr, st_devs=1):
    """

    Parameters
    ----------
    arr : np.array
        numpy array on which enhancements need to be identified.
    st_devs : float, optional
        The amount of standard deviations a peak value has to differ from the bacground concentration. The default is 1.

    Returns
    -------
    Masked array (int32), the size of input array with enhancement (1) or not (0).

    """
        
    # Step 1: Note the average and standard deviations for whole input array
    arr[arr == 0] = 'nan' # Change all zeros to 'nan' so they won't be taken into account
    average = np.nanmean(arr)
    stdev = st_devs*np.nanstd(arr) 
    arr = np.nan_to_num(arr) # Change 'nan's' back to zeros
    
    # Step 2: Isolate values that are above average
    arr[arr < average] = 0
    #average_filtered_background = np.nanmean(filtered_background)
    
    # Step 3: Check with isolated values, if they are x st.devs higher than average (background)
    arr[arr < average+stdev] = 0
    

    arr = np.nan_to_num(arr) # Rechange all 'nan' values to 0
    arr[arr > 0] = 1 # Change all enhancements to 1 (masking them)
    return arr