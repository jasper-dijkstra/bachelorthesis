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

def CheckSurroundings(arr):
    """
    Checks for amount of direct neighbors (3x3 square) around each grid cell 

    Parameters
    ----------
    arr : input array on which checks will be performed

    
    Returns
    -------
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
    
    return neighbor




def DrawBuffer(arr, buffersize):
    """
    Function to draw buffer around all grid cells with value > 0
    
    Parameters
    ----------
    arr : input array with values to draw buffers around (values > 0)
    buffersize : amount of grid cells (~7x7km) the buffer should be

    Returns
    -------
    buffer : array of size input array, with buffers drawn around all values > 0.
        buffer value == 1

    """
    
    # First create array the size of input, to draw buffers in
    buffer = np.zeros(arr.shape)
    
    # Get indices from where buffers need to be drawn
    y_indices, x_indices = np.where(arr > 0)
    
    # Loop over all x,y indices
    for arr_index in range(len(y_indices)):
        iy = y_indices[arr_index]
        ix = x_indices[arr_index]
        
        x = 0 # start x-direction
        y = buffersize # start y-direction
        diameter = 3 - 2 * buffersize # 3 (outermost layers and center point) - diameter
        
        # keep drawing new circles until buffersize is reached:
        while (x<=y): 
            for i in range(0,x + 1):
                # try-except clause to ensure working at edges of array
                try:    
                    buffer[ix+i,iy+y] = 1 #1st octant
                    buffer[ix-i,iy+y] = 1 #2nd octant
                    buffer[ix+i,iy-y] = 1 #3rd octant
                    buffer[ix-i,iy-y] = 1 #4th octant
                    buffer[ix+x,iy+i] = 1 #1st octant
                    buffer[ix-x,iy+i] = 1 #2nd octant
                    buffer[ix+x,iy-i] = 1 #3rd octant
                    buffer[ix-x,iy-i] = 1 #4th octant
                    buffer[ix+i,iy+x] = 1 #5th octant
                    buffer[ix-i,iy+x] = 1 #6th octant
                    buffer[ix+i,iy-x] = 1 #7th octant
                    buffer[ix-i,iy-x] = 1 #8th octant
                    buffer[ix+y,iy+i] = 1 #5th octant
                    buffer[ix-y,iy+i] = 1 #6th octant
                    buffer[ix+y,iy-i] = 1 #7th octant
                    buffer[ix-y,iy-i] = 1 #8th octant
                except IndexError:
                    continue # Stop this loop when edge of map is reached
            
            # Update diameter to draw new circle
            if (diameter < 0):
                diameter = diameter + 4 * x + 6
            else:
                diameter = diameter + 4 * (x-y) + 10
                y=y-1
            x=x+1
    
    return buffer
    