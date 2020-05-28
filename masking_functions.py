# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:30:07 2020

@author: Jasper Dijkstra

This script contains functions that return as masks over TROPOMI CO data (or any 2d np.array)

Functions:
    1. land/sea mask
    2. Plume mask
    
"""

import numpy as np
import warnings
from global_land_mask import globe



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


def identify_enhancements(arr, st_devs=1):
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
    
    warnings.filterwarnings('error')
    
    
    try:
        # Step 1: Note the average and standard deviations for whole input array
        arr[arr == 0] = np.nan # Change all zeros to 'nan' so they won't be taken into account
        average = np.nanmean(arr)
        stdev = st_devs*np.nanstd(arr) 
        arr = np.nan_to_num(arr) # Change 'nan's' back to zeros
        
         # Step 2: Isolate values that are above average
        arr[arr < average] = 0
        
        # Step 3: Check with isolated values, if they are x st.devs higher than average (background)
        arr[arr < average+stdev] = 0
        
    except Warning: # Is encountered above oceans, as there is no valid data
        pass
    
    arr = np.nan_to_num(arr) # Rechange all 'nan' values to 0
    arr[arr > 0] = 1 # Change all enhancements to 1 (masking them)
    
    return arr

