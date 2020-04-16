# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:30:07 2020

@author: Jasper Dijkstra

This script contains functions to detect plumes over TROPOMI CO data (or any 2d np.array), using a moving window.

Functions:
    1. land/sea mask
    2. Carbon mask
    3. Plume mask (in the future multiple!)
    
"""

import numpy as np
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



def identify_enhancements_with_quantile(frame_array, q):
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



