# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:12:14 2020

@author: Jasper Dijkstra + 2 functions from Bram Maasakkers

"""

import os
from datetime import datetime
import calendar
import numpy as np


def DefineAndCreateDirectory(targetDirectory):
    # Check if directory exists, and else create it
    if not os.path.isdir(targetDirectory):
        os.makedirs(targetDirectory)
        
        # Add notification to log file
    
    # Make sure path ends with separator (//)
    if not targetDirectory.endswith(os.path.sep):
        targetDirectory += os.path.sep
    
    return targetDirectory


def ListCSVFilesInDirectory(inputDirectory, maxfiles=None):
    # Check if directory exists, and else create it
    if not os.path.isdir(inputDirectory):
        print('This directory does not exist!')
        # Add notification to logfile
        return

    files = []
    for file in os.listdir(inputDirectory):
        if maxfiles != None:
            if len(files) == maxfiles: break # Limit the amount of input days
        if file.endswith(".csv"): files.append(os.path.join(inputDirectory, file))
        else:
            continue
    
    return files

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
    
    Convert UTC (year, month, day, hour, minute, second[, millisecond]) to modified Julian date 2000
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



def ReclassArray(daily_data_dict, arraytype, lon_resolution = 7, lat_resolution = 7):
    """
    

    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max, day, month, year,
        and np.ndarray to be reclassified.
    arraytype : string
        key of ndarray in daily_data_dict containing values to be reclassified.
    lon_resolution : integer, optional
        grid cell target resolution (km) in longitudinal direction. The default is 7.
    lat_resolution : integer, optional
        grid cell target resolution (km) in longitudinal direction. The default is 7.
        
    NOTE! lon- and lat_resolution might be distorted the further away from the equator.
    Therefore lon- and lat_resolution in km might not be exactly the input
    
    Returns
    -------
    data_reclassed : ndarray
        ndarray with resampled values of arraytype.

    """
    # Defining data_array
    data = daily_data_dict[arraytype]
    
    # Defining boundaries
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Setting the original resolution
    lon = np.linspace(lon_min, lon_max, len(data[0]))
    lat = np.linspace(lat_min, lat_max, len(data))
    
    # Setting target resolution
    nlon_t = int(abs((lon_max-lon_min)/(lon_resolution/110)))
    nlat_t = int(abs((lat_max-lat_min)/(lat_resolution/110)))
    
    # Generate target coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    #lon, lat = np.meshgrid(lon_t, lat_t)
    
    lat_list = list()
    for iobs in range(len(lat_t)):
        
        #Calculate target pixel for the observation iobs
        ilon = ((lon_t - lon[0]) / (lon[1] - lon[0])).astype('int')
        ilat = ((lat_t[iobs] - lat[0])/(lat[1]-lat[0])).astype('int')
        
        lat_row = data[ilat,ilon]
        lat_list.append(lat_row)
     
    data_reclassed = np.array(lat_list)
    
    return data_reclassed
