# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:30:44 2020

@author: Jasper Dijkstra

This script contains a 2d moving window function over a np.ndarray

"""

import numpy as np
import masking_functions as mask


def MovingWindow(arr, window = (100,100), step = 20, treshold = 0.5):
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

    Returns
    -------
    np.array of the same size as input arrays, but masked enhancements.

    """
    
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
            #frame = mask.identify_enhancements_1(frame, q=0.9)
            #frame = mask.identify_enhancements_2(frame, q=0.9)
            frame = mask.identify_enhancements_3(frame, st_devs = 2)
            
            #Append the moving window to the count and field arrays:
            field[y:y+window[1],x:x+window[0]] = frame[0:window[1],0:window[0]] + field[y:y+window[1],x:x+window[0]]
            count[y:y+window[1],x:x+window[0]] = count_frame[0:window[1],0:window[0]] + count[y:y+window[1],x:x+window[0]]
    
            y += step
        x += step
    
    # Decide upon real plumes
    field = field/count
    field[field >= treshold] = 1.
    field[field < treshold] = 0.
    
    return field

def CheckForSurroundings(arr, neighbors=1):
    """
    Checks for amount of direct neighbors (3x3 square) around each grid cell 

    Parameters
    ----------
    arr : input array on which checks will be performed
    neighbors : int between 0-9 (default = 1)
        minimum amount of neighbors required to be considered a plume
    
    Returns
    -------
    arr : input array without grid cells with too few neighbors
    neighbor : array of size arr, with the amount of neighbors of each grid cell

    """
    # Create array of size input array
    neighbor = np.zeros(arr.shape, dtype=float)
    
    # Check for neighbours in {}-direction:
    neighbor[:-1] += arr[1:]        # North
    neighbor[:-1,1:] += arr[1:,:-1] # North East
    neighbor[:,1:] += arr[:,:-1]    # East
    neighbor[1:,1:] += arr[:-1,:-1] # South East
    neighbor[1:] += arr[:-1]        # South
    neighbor[1:,:-1] += arr[:-1,1:] # South West
    neighbor[:,:-1] += arr[:,1:]    # West
    neighbor[:-1,:-1] += arr[1:,1:] # North West
    
    # Removing nuisances by making sure there is at least one neighbour
    arr[neighbor < neighbors] = 0
        
    return {'array' : arr, 'neighbors' : neighbor}