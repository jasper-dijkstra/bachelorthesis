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

from datetime import datetime
import numpy as np
import pandas as pd
import csv

reference_date = datetime(2010, 1, 1, 0, 0, 0)

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

def reading_csv_as_df(csv_file):
    """
    --------
    Parameters:
        csvfile: path to .csv file containing monthly TROPOMI data on CO
    
    Returns:
        list with all grid cells that fall within bbox for input within timeframe of csv file
    """
    
    CO_header = ['lat0', 'lat1', 'lat2', 'lat3', 'lat', 'lon0', 'lon1', 'lon2', 'lon3', 'lon', 'xco', 'xco_unc', 'xco_ppb', 'xco_ppb_unc', 'qa', 'weekday', 'day', 'month', 'aerosol_opthick', 'aerosol_layer', 'orbitnr', 'time']
    COdata = pd.read_csv(csv_file, header=None, names=CO_header)
    return COdata


def reading_csv_as_nparray(csvfile, bbox, target_lon, target_lat):
    """
    
    Parameters
    ----------
    csvfile : string
        Path to input csv file.
    bbox : list
        List containing the minimum and maximum longitudes (x) and latitudes (y) for area of interest [lat_min, lat_max, lon_min, lon_max]
    
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
    
    COdata = reading_csv_as_df(csvfile) # Reading Australia csv, as pd.dataframe
    
    # Retrieving day, month and year of dataset
    day = COdata.at[0, 'day']
    month = COdata.at[0, 'month']
    seconds = reference_date.timestamp() + int(COdata.at[0, 'time'])
    year = datetime.fromtimestamp(seconds).year

        
    target = 'xco_ppb' # Set target for map creation
    
    # Target Resolution
    nlon_t=target_lon
    nlat_t=target_lat
    
    # Initiate arrays that save the observation totals for every pixel
    field_t = np.zeros((nlat_t,nlon_t))
    count_t = np.zeros((nlat_t,nlon_t))
    
    delta_lon = int(lon_max - lon_min) # Calculate Target Longitudal Range
    delta_lat = int(lat_max - lat_min) # Calculate Target Latitudal Range
    
    # reduce data resolution to target resolution
    for iobs in range(len(COdata)):
        
        #Find observation coordinates
        lon_obs = COdata.at[iobs, 'lon']
        lat_obs = COdata.at[iobs, 'lat']
                
        #Calculate target pixel for the observation iobs
        ilon = np.int((lon_obs - lon_min)* nlon_t / delta_lon)
        ilat = np.int((lat_obs - lat_min)* nlat_t / delta_lat)
        
        if COdata.at[iobs, target] > 2.5e3: continue
        field_t[ilat,ilon] += COdata.at[iobs, target]
        count_t[ilat,ilon] += 1
    
    idx = (count_t > 0)
    field_t[idx] = field_t[idx]/count_t[idx]
    returndict = {'day':day, 'month':month, 'year':year, 'lon_min':lon_min, 'lon_max':lon_max, 'lat_min':lat_min, 'lat_max':lat_max, 'target_lon':nlon_t, 'target_lat':nlat_t, 'count_t':count_t, 'CO_ppb':field_t}
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



