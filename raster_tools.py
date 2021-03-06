# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:18:23 2020

@author: Jasper Dijkstra

Raster functions, to be applied on 2D np.arrays

All functions in this file:
    - MovingWindow() - Apply moving window over 2D array
    - CountNeighbors() - Defines the amount of direct neighbors (3x3 square) around each grid cell 
    - DrawCircularBuffer() - Apply a midpoint circle algorithm
    - ResampleArray() - Resample array to a different resolution.

"""

import numpy as np



def MovingWindow(arr, function, window = (100,100), step = 20, st_devs=1):
    """
    Function to apply moving window on 2D np.array <arr>. 
    <function> will be applied on <window>, that moves with steps of size <step>
    
    
    Parameters
    ----------
    arr : np.array
        2D np.array over which a moving window with function <function> will be applied.
    function : function
        Function that takes an array as input. 
        Function will be calculated for each window.
    window : tuple (x,y) with integers.
        Tuple that describes the size of the moving window (x, y). 
        Default is (100,100)
    step : Integer.
        The amount of pixels the window has to step. Default is 20

    Returns
    -------
    np.array of the same size as input array, but <function> applied to each <window>.

    """
    
    # Assert window size is larger than or equal to step size, else values would be skipped
    assert window[0] >= step and window[1] >= step, 'step size cannot be larger than window size!'
    
    # Define the np.arrays where the data has to be appended on
    field = np.zeros((len(arr),len(arr[0])))
    count = np.zeros((len(arr),len(arr[0])))
    
    # Define frame average, list with the average CO concentration per frame
    average = list()
    
    x = 0
    for xdir in range(int(len(arr[0])/step)):
        y = 0
        for ydir in range(int(len(arr)/step)):
            frame = arr[y:y+window[1],x:x+window[0]] # Moving frame 
            count_frame = np.full((len(frame),len(frame[0])), 1) # Count frame
            
            # Analysis applied onto frame
            frame, frame_average = function(frame, st_devs)

            #Append the moving window to the count and field arrays:
            field[y:y+window[1],x:x+window[0]] = frame[0:window[1],0:window[0]] + field[y:y+window[1],x:x+window[0]]
            count[y:y+window[1],x:x+window[0]] = count_frame[0:window[1],0:window[0]] + count[y:y+window[1],x:x+window[0]]
            
            # Append the frame average to list
            average.append(frame_average)
            
            y += step
        x += step
    
# =============================================================================
#     try:
#         field = field/count#.astype(int)
#     except RuntimeWarning:
# =============================================================================
        
    field = np.divide(field, count, where=(count != 0))
    
    # Return the average of the average list
    avg = np.nanmean(np.array(average))
    
    return field, avg


def CountNeighbors(arr):
    """
    Defines the amount of direct neighbors (3x3 square) around each grid cell 

    Parameters
    ----------
    arr : input array on which neighbors will be performed
    
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


def DrawCircularBuffer(arr, radius):
    """
    Function to apply a midpoint circle algorithm on all grid cells of <arr> with value > 0
    
    Parameters
    ----------
    arr : input array with values to draw buffers around (values > 0)
    radius : radius (in grid cells) the buffer should be

    Returns
    -------
    buffer : array of size input array, with buffers drawn around all values > 0.
        buffer value == 1

    """
    
    # Function used to append rastervalues to correct index (if exists)
    def appendone(buffer,y,x):
        try:
            buffer[y,x] = 1
        except IndexError:
            pass
        return
    
    # First create array the size of input, to draw buffers in
    buffer = np.zeros(arr.shape)
    
    # Get indices from where buffers need to be drawn
    x_indices, y_indices = np.where(arr > 0)
    
    # Loop over all x,y indices
    for arr_index in range(len(y_indices)):
        iy = y_indices[arr_index]
        ix = x_indices[arr_index]
        
        x = 0 # start x-direction
        y = radius # start y-direction
        diameter = 3 - 2 * radius # 3 (outermost layers and center point) - diameter
        
        # keep drawing new circles until buffersize is reached:
        while (x<=y): 
            for i in range(0,x + 1):
                # Append the value of one in seperate try-except function,
                # to avoid problems at map edges
                appendone(buffer,ix+i,iy+y) #1st octant
                appendone(buffer,ix-i,iy+y) #2nd octant
                appendone(buffer,ix+i,iy-y) #3rd octant
                appendone(buffer,ix-i,iy-y) #4th octant
                appendone(buffer,ix+x,iy+i) #1st octant
                appendone(buffer,ix-x,iy+i) #2nd octant
                appendone(buffer,ix+x,iy-i) #3rd octant
                appendone(buffer,ix-x,iy-i) #4th octant
                appendone(buffer,ix+i,iy+x) #5th octant
                appendone(buffer,ix-i,iy+x) #6th octant
                appendone(buffer,ix+i,iy-x) #7th octant
                appendone(buffer,ix-i,iy-x) #8th octant
                appendone(buffer,ix+y,iy+i) #5th octant
                appendone(buffer,ix-y,iy+i) #6th octant
                appendone(buffer,ix+y,iy-i) #7th octant
                appendone(buffer,ix-y,iy-i) #8th octant

            # Update diameter to draw new circle
            if (diameter < 0):
                diameter = diameter + 4 * x + 6
            else:
                diameter = diameter + 4 * (x-y) + 10
                y=y-1
            x=x+1
    
    return buffer


def ResampleArray(bbox, array, lon_resolution, lat_resolution):
    """
    Resampling <array> to a different resolution.

    Parameters
    ----------
    bbox : list
        [lat_min, lat_max, lon_min, lon_max] of output extents
    array : np.array
        np.ndarray containing values to be reclassified.
    lon_resolution : integer
        grid cell target resolution (km) in longitudinal direction.
    lat_resolution : integer
        grid cell target resolution (km) in longitudinal direction.
        
    NOTE! lon- and lat_resolution might be distorted the further away from the equator.
    Therefore lon- and lat_resolution in km might not be exactly the input
    
    Returns
    -------
    data_reclassed : ndarray
        ndarray with resampled values of array.

    """
    # Defining data_array
    data = array
    
    # Defining boundaries
    lat_min = bbox[0]
    lat_max = bbox[1]
    lon_min = bbox[2]
    lon_max = bbox[3]
    
    # Setting the original resolution
    lon = np.linspace(lon_min, lon_max, len(data[0]))
    lat = np.linspace(lat_min, lat_max, len(data))
    
    # Setting target resolution
    nlon_t = int(abs((lon_max-lon_min)/(lon_resolution/110)))
    nlat_t = int(abs((lat_max-lat_min)/(lat_resolution/110)))
    
    # Generate target coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    
    lat_list = list()
    for iobs in range(len(lat_t)):
        
        #Calculate target pixel for the observation iobs
        ilon = ((lon_t - lon[0]) / (lon[1] - lon[0])).astype('int')
        ilat = ((lat_t[iobs] - lat[0])/(lat[1]-lat[0])).astype('int')
        
        lat_row = data[ilat,ilon]
        lat_list.append(lat_row)
     
    data_reclassed = np.array(lat_list)
    
    return data_reclassed


