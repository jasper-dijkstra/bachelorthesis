# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:12:14 2020

@author: Jasper Dijkstra + 2 functions from Bram Maasakkers

This script contains functions to:
    1. Check if directory exists, if not create it
    2. List all files in directory
    3. Export data to csv file
    4. Get current time (year, month, day, hour, minute, second)
    5. Convert UTC to modified Julian date 2010 ((C) Bram Maasakkers)
    6. Convert modified Julian date 2010 to UTC ((C) Bram Maasakkers)

"""

import os
import csv
from datetime import datetime
import calendar
import numpy as np

# ==============================================
# OS FUNCTIONS
# ==============================================

def DefineAndCreateDirectory(targetDirectory):
    """
    Check if directory exists, and else create it 
    """
    
    if not os.path.isdir(targetDirectory):
        os.makedirs(targetDirectory)
    
    # Make sure path ends with separator (//)
    if not targetDirectory.endswith(os.path.sep):
        targetDirectory += os.path.sep
    
    return targetDirectory



def ListFilesInDirectory(inputDirectory, extension='.csv', maxfiles=None):
    """
    list all files of given extension (default is '.csv') in a directory.
    """
    
    # Check if directory exists, else raise AssertionError
    assert os.path.isdir(inputDirectory), 'This directory does not exist!'

    files = []
    for file in os.listdir(inputDirectory):
        if maxfiles != None:
            if len(files) == maxfiles: break # Limit the amount of input days
        if file.endswith(extension): files.append(os.path.join(inputDirectory, file))
        else:
            continue
    
    return files


def ExportAsCSV(csv_out_path, data):
    """
    export data (lists) as a csv file

    Parameters
    ----------
    csv_out_path : string
        Path to output csv file.
    data : list
        list that contains lists with data to be used as output.

    Returns
    -------
    None.

    """
    with open(csv_out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(data)
        f.close()
    return




# ==============================================
# TIME CONVERSION FUNCTIONS
# ==============================================

def GetCurrentTime():
    now = datetime.now()
    year = now.strftime("%Y")
    month = now.strftime("%m")
    day = now.strftime("%d")
    hour = now.strftime("%H")
    minute = now.strftime("%M")
    second = now.strftime("%S")
    
    return {'year' : year, 'month' : month, 'day' : day, 'hour' : hour, 'minute' : minute, 'second' : second}


def UTCtoModifiedJulianDate(year, month, day, hour, minute, second, millisecond=0):
    """
    FUNCTION BY BRAM MAASAKKERS (C)
    
    Convert UTC (year, month, day, hour, minute, second[, millisecond]) to modified Julian date 2010
    This function is vector-safe.
    
    Parameters:
        year: the year
        month: the month
        day: the day
        hour: the hour
        minute: the minute
        second: the second
        millisecond: the millisecond (default: 0)
    Return value:
        result: the time/date converted to modified Julian date 2010
    """

    t2000 = calendar.timegm((2010, 1, 1, 0, 0, 0)) # seconds in epoch 1970 corresponding to 1 Jan 2000
    if isinstance(year, np.ndarray): # if input is vectors
        is_scalar = False
        if isinstance(millisecond,int): # if millisecond is an int, which means that the optional argument has presumably not been set
            millisecond = np.zeros(year.shape, dtype=np.int) # set millisecond to an array of matching size
        elif any(millisecond != 0) and any(second%1 != 0): # if both millisecond and fractional second are given
            print("Warning: both milliseconds and fractional seconds given! Ignoring fractional seconds.")
            second = np.floor(second).astype(np.int)
    else: # input is scalars
        is_scalar = True
        if millisecond != 0 and second%1 != 0: # if both millisecond and fractional second are given
            print("Warning: both milliseconds and fractional seconds given! Ignoring fractional seconds.")
            second = np.int(np.floor(second)) # cut off fractional seconds
        year = np.array([year]) # convert to arrays
        month = np.array([month])
        day = np.array([day])
        hour = np.array([hour])
        minute = np.array([minute])
        second = np.array([second])
        millisecond = np.array([millisecond])
    result = np.zeros(year.shape, dtype=np.float64) # initialise field for result
    for ind in range(len(year)): # loop over entries of vectors (one entry for scalars)
        t = calendar.timegm((year[ind], month[ind], day[ind], hour[ind], minute[ind], second[ind])) # convert to seconds in epoch 1970
        result[ind] = ( np.float64(t - t2000) + np.float64(millisecond[ind])/1000.0 ) / 86400.0 # convert to (fractional) days since 1 Jan 2000
    if is_scalar: # if input has been scalars
        result = result.item() # convert result back to scalar
    return result


def ModifiedJulianDatetoUTC(mjd):
    """
    FUNCTION BY BRAM MAASAKKERS (C)
    
    Convert modified Julian date 2010 to UTC
    This function is vector-safe.
    
    Parameters:
        mjd: time/date as modified Julian date 2000
    Return value:
        result: a dictionary with the entries "year", "month", "day", "hour", "minute", "second", "millisecond", "fractional_year"
    """
    import calendar
    import time
    
    global t
    global gmt
    global t2010
    
    t2010 = calendar.timegm((2010, 1, 1, 0, 0, 0)) # seconds in epoch 1970 corresponding to 1 Jan 2010
    
    # Vectorize inputs
    if isinstance(mjd, (np.ndarray,list)): # if input is a vector
        is_scalar = False
        if isinstance(mjd,list):
            mjd = np.array(mjd)
    else: # input is a scalar
        is_scalar = True
        mjd = np.array([mjd])  # convert to array
        
    t = mjd + t2010 # compute seconds since epoch 1970
    gmt = np.zeros((len(mjd),9), dtype=np.int) # initialise field for intermediate result
    fractional_year = np.zeros(len(mjd), dtype=np.double)
    day_of_year = np.zeros(len(mjd), dtype=np.double)
    for ind in range(len(t)): # loop over entries of vector (may be one entry)
        gmt[ind,:] = np.array(time.gmtime(t[ind]))
        mjd_year_begin = UTCtoModifiedJulianDate(gmt[ind,0], 1, 1, 0, 0, 0)
        day_of_year[ind] = (mjd[ind] - mjd_year_begin)
        fractional_year[ind] = float(gmt[ind,0]) + day_of_year[ind] / (366.0 if calendar.isleap(gmt[ind,0]) else 365.0)
    if is_scalar:
        result={"year":gmt[0,0],"month":gmt[0,1], "day":gmt[0,2], "hour":gmt[0,3], "minute":gmt[0,4], "second":gmt[0,5], "millisecond":np.int(np.round(t[0]%1*1000.0)), "fractional_year":fractional_year[0], "day_of_year":day_of_year[0]}
    else:
        result={"year":gmt[:,0],"month":gmt[:,1], "day":gmt[:,2], "hour":gmt[:,3], "minute":gmt[:,4], "second":gmt[:,5], "millisecond":np.round(np.array(t%1)*1000.0).astype(np.int), "fractional_year":fractional_year, "day_of_year":day_of_year}
    return result





    