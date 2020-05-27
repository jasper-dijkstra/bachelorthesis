# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:25:51 2020

@author: Jasper Dijkstra, edited from scripts from M. Riet, S. Houweling and I. Dekker

This script contains standalone functions to:
    1. filter TROPOMI csv CO data by day
    2. filter TROPOMI csv CO data by area (bbox)
    3. reading csv as pandas dataframe
    4. reading csv as np.array
    5. export to csv file


"""

import numpy as np
import pandas as pd
import csv

import utilities as ut

#----------------------------------
# FILTERING DATA FUNCTIONS
#----------------------------------

def filter_csv_by_day(csvfile, day):
    """
    --------
    Parameters:
        csvfile: path to .csv file containing monthly TROPOMI data on CO
        day: Day of the data you want to filter
    
    Returns:
        list with all grid cells that fall within selected day
    """
    
    output_csv = []
    with open(csvfile, newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if int(row[16]) == int(day):
                output_csv.append(row) # Appending all output data
        f.close()
    return output_csv



def filter_csv_by_bbox(csvfile, bbox):
    """
    --------
    Parameters:
        csvfile: path to .csv file containing monthly TROPOMI data on CO
        bbox: List containing the minimum and maximum longitudes (x) and latitudes (y) for area of interest [lat_min, lat_max, lon_min, lon_max]
    
    Returns:
        list with all grid cells that fall within bbox for input within timeframe of csv file
    """
    
    # Defining boundaries
    latmin = bbox[0]
    latmax = bbox[1]
    lonmin = bbox[2]
    lonmax = bbox[3]
    
    output_csv = []
    with open(csvfile, newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            y = float(row[4]) # y (latitude) coordinate of center gridcell
            x = float(row[9]) # x (longitude) coordinate of center gridcell
            if x > lonmin and x < lonmax and y > latmin and y < latmax:
                output_csv.append(row) # Appending all output data
        f.close()
    return output_csv


#----------------------------------
# IMPORTING AND EXPORTING CSV FILES
#----------------------------------

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
    qa = np.zeros((target_lat, target_lon))
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
            qa[ilat,ilon] += qa_1d[iobs]
            hour[ilat, ilon] += hour_1d[iobs]
            count_t[ilat,ilon] += 1
        else:
            pass
    

    # Removing data with a too high uncertainty
    #xco_ppb[qa < max_unc] = 0
    #hour[qa < max_unc] = 0
    #count_t[qa < max_unc] = 0
    
    # Apply correction for cells that overlap
    idx = (count_t > 0)
    xco_ppb[idx] = xco_ppb[idx]/count_t[idx]
    #qa[idx] = qa[idx]/count_t[idx]
    hour[idx] = hour[idx]/count_t[idx]
    
    # Change all 'no data' values (count_t == 0) to np.nan
    hour[count_t == 0] = np.nan

   
    returndict = {'day':day, 'month':month, 'year':year, \
                  'lon_min':lon_min, 'lon_max':lon_max, 'lat_min':lat_min,\
                      'lat_max':lat_max, 'target_lon':target_lon, 'target_lat':target_lat,\
                          'count_t':count_t, 'CO_ppb':xco_ppb, 'hour': hour}

    return returndict


def export_as_csv(csv_out_path, data):
    """

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


