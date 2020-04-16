# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:30:44 2020

@author: jaspd
"""

import numpy as np
import masking_functions as mask
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
            frame = mask.identify_enhancements_with_quantile(frame, q)
            
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