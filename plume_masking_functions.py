# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:30:07 2020

@author: Jasper Dijkstra

This script contains functions to detect plumes over TROPOMI CO data (or any 2d np.array), using a moving window.

Functions:
    1. moving window
    2. Plume mask
    
"""

import numpy as np
from datetime import datetime

def moving_window(arr, window = (100,100), step = 20, treshold = 0.5, q = 0.9):
    """
    
    Parameters
    ----------
    arr : np.array
        2D np.array with TROPOMI data on CO total column concentration (ppb).
    window : tuple (x,y) with integers.
        Tuple that describes the size of the moving window (x, y). Default is (100,100)
    step : Integer.
        The amount of pixels the window has to step. Default is 20
    treshold : float (0-1).
        percentage of minimum amount of windows a peak has to be detected in. Default is 0.5
    q : float (0-1)
        The q-th quantile a value has to be, to be considered a peak. Default is 0.9

    Returns
    -------
    np.array of the same size as input arrays, but masked enhancements.

    """
    start = datetime.now()
    
    # Define the np.arrays where the data has to be appended on
    field = np.zeros((len(arr),len(arr[0])))
    count = np.zeros((len(arr),len(arr[0])))
    
    x = 0
    for xdir in range(int(len(arr[0])/step)):
        y = 0
        for ydir in range(int(len(arr)/step)):
            frame = arr[y:y+window[1],x:x+window[0]] # Moving frame 
            count_frame = np.full((len(frame),len(frame[0])), 1) # Count frame
            
            # Analysis applied onto frame
            frame = identify_enhancements_with_quantile(frame, q)
            
            #Append the moving window to the count and field arrays:
            field[y:y+window[1],x:x+window[0]] = frame[0:window[1],0:window[0]] + field[y:y+window[1],x:x+window[0]]
            count[y:y+window[1],x:x+window[0]] = count_frame[0:window[1],0:window[0]] + count[y:y+window[1],x:x+window[0]]
    
            y += step
        x += step
    
    # Decide upon real plumes
    field = field/count
    field[field >= treshold] = 1.
    field[field < treshold] = 0.
    
    print('Total elapsed time for moving window: {}'.format(datetime.now()-start))
    
    return field

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

