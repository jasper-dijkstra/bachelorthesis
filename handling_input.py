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
    6. read data from Global Fire Emissions Database (hdf5)

"""

import os
import numpy as np
import pandas as pd
import csv
import h5py

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


def reading_csv_as_nparray(csvfile, bbox, target_lon, target_lat, max_unc):
    """
    
    Parameters
    ----------
    csvfile : string
        Path to input csv file.
    bbox : list
        List containing the minimum and maximum longitudes (x) and latitudes (y) for area of interest [lat_min, lat_max, lon_min, lon_max]
    target_lon
    target_lat
    max_unc : maximum uncertainty in TROPOMI measurements
    
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
    # Create dict form Julian Date to UTC
    time_dict = ut.ModifiedJulianDatetoUTC(int(COdata.at[0, 'time']))
    #hour = time_dict['hour'] # This can be an enhancement for meteodata!
    year = time_dict['year']
        
    co_ppb = 'xco_ppb'      # Set target for map creation
    uncertainty = 'qa' #'xco_ppb_unc'
    
    # Target Resolution
    nlon_t=target_lon
    nlat_t=target_lat
    
    # Initiate arrays that save the observation totals for every pixel
    field_t = np.zeros((nlat_t,nlon_t))
    count_t = np.zeros((nlat_t,nlon_t))
    unc_t = np.zeros((nlat_t,nlon_t))
    
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
        
        if COdata.at[iobs, co_ppb] > 2.5e3: continue
        field_t[ilat,ilon] += COdata.at[iobs, co_ppb]
        unc_t[ilat,ilon] += COdata.at[iobs, uncertainty]
        count_t[ilat,ilon] += 1
    
    idx = (count_t > 0)
    field_t[idx] = field_t[idx]/count_t[idx]
    
    # Removing data with a too high uncertainty
    field_t[unc_t > max_unc] = 0
    count_t[unc_t > max_unc] = 0
    
    returndict = {'day':day, 'month':month, 'year':year, 'uncertainty':unc_t,\
                  'lon_min':lon_min, 'lon_max':lon_max, 'lat_min':lat_min,\
                      'lat_max':lat_max, 'target_lon':nlon_t, 'target_lat':nlat_t,\
                          'count_t':count_t, 'CO_ppb':field_t}#, 'hour': hour}
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


def ReadGFED(daily_data_dict):
    """
    This function reads HDF5 data from the Global Fire Emissions Database (GFED),
    and then resample it to the desired resolution. GFED data from 2017 onwards,
    is only available as Beta version. Therefore this function makes no distinction
    in plume intensity.
    
    Parameters
    ----------
    daily_data_dict : dictionary
        dictionary containing daily_data per day [<day>].
    
    Returns
    -------
    array of same resoultion and scope as described in daily_data_dict, with GFED values
    
    """
    
    # STEP 1: OBTAIN DATA FROM DAILY_DATA_DICT
    # Define lat/lon ranges of target
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Define timescope
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    
    
    # STEP 2: READ THE DATA FROM THE HDF5 FILE
    # Get Path of current script
    full_path = str(os.path.realpath('__file__'))
    path, _ = os.path.split(full_path)
    
    # Navigate to the correct hdf5 folder
    gfed_filename = str(path + r'\GFED_fire_emissions_mask\GFED4.1s_{}_beta.hdf5'.format(year))
    f = h5py.File(gfed_filename, mode='r')
    
    # Get daily fire and lat/lon data from hdf5 file
    monthly_data = f['/emissions/{}/C'.format(month)][:]
    daily_fraction = f['/emissions/{}/daily_fraction/day_{}'.format(month, day)][:]
    data = np.multiply(monthly_data, daily_fraction)
    lat_hdf = f['/lat'][:]
    lat_hdf = lat_hdf[:, 0]
    lon_hdf = f['/lon'][:]
    lon_hdf = lon_hdf[0]
    
    
    # STEP 3: RESAMPLE TO TARGET RESOLUTION
    # Defining the lat and longitudinal ranges of output
    latrange = np.linspace(lat_min, lat_max, len(daily_data_dict['CO_ppb']))
    lonrange = np.linspace(lon_min, lon_max, len(daily_data_dict['CO_ppb'][0]))
    
    gfed = list()
    for iobs in range(len(latrange)):
                
        #Calculate target pixel for the observation iobs
        ilon = ((lonrange - lon_hdf[0]) / (lon_hdf[1] - lon_hdf[0])).astype('int')
        ilat = ((latrange[iobs] - lat_hdf[0])/(lat_hdf[1]-lat_hdf[0])).astype('int')
        
        # Append data to correct lat/lon indices
        gfed_row = data[ilat,ilon]
        gfed.append(gfed_row)
    gfed = np.array(gfed)
    
    # Set all values > 0 to 1
    #gfed[gfed > 0] = 1
    
    return gfed
