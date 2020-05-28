# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:25:51 2020

@author: Jasper Dijkstra, edited from scripts from M. Riet, S. Houweling and I. Dekker


This script contains functions to handle daily TROPOMI csv files.
- Data required for plume detection algorithm is extracted and returned


"""

import numpy as np
import pandas as pd

# Local imports
import utilities as ut


def CSVtoDf(csv_file, bbox):
    """
    --------
    Parameters:
        csvfile: path to .csv file containing monthly TROPOMI data on CO
    
    Returns:
        list with all grid cells that fall within bbox for input within timeframe of csv file
    """
    
    # Dataframe header values
    CO_header = ['lat0', 'lat1', 'lat2', 'lat3', 'lat', 'lon0', 'lon1', 'lon2', 'lon3', 'lon', 'xco', 'xco_unc', 'xco_ppb', 'xco_ppb_unc', 'qa', 'weekday', 'day', 'month', 'aerosol_opthick', 'aerosol_layer', 'orbitnr', 'time']
    
    # Read the csv file as Pandas DataFrame
    COdata = pd.read_csv(csv_file, header=None, names=CO_header)
    
    # Identify if observations falls within latlon range
    COdata['valid'] = (COdata['lat'] >= bbox[0]) & (COdata['lat'] < bbox[1]) \
        & (COdata['lon'] >= bbox[2]) & (COdata['lon'] < bbox[3])
    
    # Remove all values that do not fall within latlon range
    COdata = COdata.drop(COdata[COdata['valid'] == False].index)
    del COdata['valid']
    
    return COdata


def CSVtoArray(csvfile, bbox, target_lon, target_lat):
    """
    
    Parameters
    ----------
    csvfile : string
        Path to input csv file.
    bbox : list
        List containing the minimum and maximum longitudes (x) and latitudes (y) for area of interest [lat_min, lat_max, lon_min, lon_max]
    target_lon
    target_lat
    
    Returns
    -------
    returnlist : list
        list containing: lon_min, lon_max, lat_min, lat_max, nlon_t, nlat_t, count_t, field_t.

    """
    
    # Defining boundaries
    lat_min = bbox[0]
    lat_max = bbox[1]
    lon_min = bbox[2]
    lon_max = bbox[3]
    
    
    COdata = CSVtoDf(csvfile, bbox) # Reading csv as pd.dataframe
    COdata = COdata.to_numpy() # Transforming pd.dataframe to np.array

    # Convert TROPOMI timestamp to UTC date, time floats
    time = ut.ModifiedJulianDatetoUTC(COdata[:,21])
    
    # Get year, month, day from time dictionary (not spatially explicit)
    year = time['year'][0]
    month = time['month'][0]
    day = time['day'][0]
    
    # Spatially explicit data from TROPOMI
    lat = COdata[:,4]
    lon = COdata[:,9]
    xco_ppb_1d = COdata[:,12]
    qa_1d = COdata[:,14]
    hour_1d = time['hour']
    
    
    # Give indices with a too high a 'nodata' value
    qa_indices = np.where((qa_1d <= 0.1) | (qa_1d >= 2.0))[0]
    xco_ppb_1d[qa_indices] = np.nan
    
    # Initiate arrays that save the observation totals for every pixel
    xco_ppb = np.zeros((target_lat, target_lon))
    hour = np.zeros((target_lat, target_lon))
    count_t = np.zeros((target_lat, target_lon))


    # reduce data resolution to target resolution
    for iobs in range(len(xco_ppb_1d)):
        if iobs != np.nan: # Make sure only valid observations are placed in grid
            #Calculate target pixel for the observation iobs    
            ilon = np.int((lon[iobs] - lon_min)* target_lon / (lon_max - lon_min))
            ilat = np.int((lat[iobs] - lat_min)* target_lat / (lat_max - lat_min))
            
            # Append iobs to correct ilat, ilon (indices) in 2d array
            xco_ppb[ilat,ilon] += xco_ppb_1d[iobs]
            hour[ilat, ilon] += hour_1d[iobs]
            count_t[ilat,ilon] += 1
        else:
            pass
    
    # Apply correction for cells that overlap
    idx = (count_t > 0)
    xco_ppb[idx] = xco_ppb[idx]/count_t[idx]
    hour[idx] = hour[idx]/count_t[idx]
    
    # Change all 'no data' values (count_t == 0) to np.nan
    hour[count_t == 0] = np.nan
   
    returndict = {'day':day, 'month':month, 'year':year, \
                  'lon_min':lon_min, 'lon_max':lon_max, 'lat_min':lat_min,\
                      'lat_max':lat_max, 'target_lon':target_lon, 'target_lat':target_lat,\
                          'count_t':count_t, 'CO_ppb':xco_ppb, 'hour': hour}

    return returndict





