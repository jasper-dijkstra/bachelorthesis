# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 21:38:48 2020

@author: Jasper Dijkstra

This script contains functions to:
    1. Download ERA5 meteodata on u- and v-winds for a desired atmospheric pressure
        level, time- and aerial scope.
    2. Process the data obtained with function 1, to np arrays the size of input
    
"""

from netCDF4 import Dataset
import cdsapi
import numpy as np
import os

# Local imports
import utilities as ut


def DownloadERA5(daily_data_dict, pressure_level=850):
    """
    This function downloads ERA5 data making use of the CDS API from Copernicus,
    as netCDF4 files for the time- and aerial scope defined in daily_data_dict
    
    Parameters
    ----------
    daily_data_dict : dictionary
        dictionary with TROPOMI data.
    pressure_level : integer, string or float
        atmospheric pressure level (hPa) of the winddata. The default is 850.

    Returns
    -------
    downloaded .nc file in: (<working_dir>/store_meteo_data).

    """
    # STEP 1: Preparing input for API request:
    days = list()
    months = list()
    years = list()
    for day in daily_data_dict:
        day_nr = str(daily_data_dict[day]['day'])
        month_nr = str(daily_data_dict[day]['month'])
        year_nr = str(daily_data_dict[day]['year'])
        days.append(day_nr)
        months.append(month_nr)
        years.append(year_nr)
    
    
    curr_directory = os.getcwd()
    out_dir = os.path.join(curr_directory, r'store_meteo_data')
    ut.DefineAndCreateDirectory(out_dir)
    filename = 'ERA5_{}{}{}-{}{}{}.nc'.format(min(months), min(days), min(years), max(months), max(days), max(years))
    
    if os.path.isfile(os.path.join(out_dir, filename)) == True:
        print('ERA5 meteodata already exists for the period: {}/{}/{}-{}/{}/{}!'.format(min(months), min(days), min(years), max(months), max(days), max(years)))
        print('Check if extents of data are also correct!')
        return os.path.join(out_dir, filename)
    
    lat_min = daily_data_dict[0]['lat_min']
    lat_max = daily_data_dict[0]['lat_max']
    lon_min = daily_data_dict[0]['lon_min']
    lon_max = daily_data_dict[0]['lon_max']
    
    # STEP 2: API request:
    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'area': [lat_min, lon_min, lat_max, lon_max],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'day': days,
            'month': months,
            'year': years,
            'pressure_level': str(pressure_level),
            'variable': ['u_component_of_wind', 'v_component_of_wind'],
        },
        os.path.join(out_dir, filename))
    
    return os.path.join(out_dir, filename)


def ProcessMeteo(daily_data_dict, path):
    """
    Function that processes data obtained via DownloadERA5() to a u- and v- wind
    array the same resolution as TROPOMI.
    
    Parameters
    ----------
    daily_data_dict: dictionary with TROPOMI data
    path: path (str) to downloaded ERA5 (.nc) file
    
    ----------
    Returns: updated daily_data_dict with u- and v-winds array

    """
    #Defining lat- and longitudinal ranges
    lat_min = daily_data_dict[0]['lat_min']
    lat_max = daily_data_dict[0]['lat_max']
    lon_min = daily_data_dict[0]['lon_min']
    lon_max = daily_data_dict[0]['lon_max']
    
    fid = Dataset(path)
    
    lon_meteo = fid.variables['longitude'][:]
    lat_meteo = fid.variables['latitude'][:]
    time = fid.variables['time'][:]
    
    index = 0
    for i in range(int(len(time)/24)):
        # Reading the daily mean of u and v
        u_wind = fid.variables['u'][index:index+24]
        u_wind = np.mean(u_wind, axis=0)
        v_wind = fid.variables['v'][index:index+24]
        v_wind = np.mean(v_wind, axis=0)
        index += 24 # 24 hours
    
        # Defining the lat and longitudinal ranges of u and v arrays
        latrange = np.linspace(lat_min, lat_max, len(daily_data_dict[0]['CO_ppb']))
        lonrange = np.linspace(lon_min, lon_max, len(daily_data_dict[0]['CO_ppb'][0]))
        
        # Reclassifying u and v wind to target resolution
        u_wind_list = list()
        v_wind_list = list()
        for iobs in range(len(latrange)):   
            #Calculate target pixel for the observation iobs
            ilon = ((lonrange - lon_meteo[0]) / (lon_meteo[1] - lon_meteo[0])).astype('int')
            ilat = ((latrange[iobs] - lat_meteo[0])/(lat_meteo[1] - lat_meteo[0])).astype('int')
            
            # Append data to correct lat/lon indices
            u_row = u_wind[ilat,ilon]
            u_wind_list.append(u_row)
            v_row = v_wind[ilat,ilon]
            v_wind_list.append(v_row)
            
        v_wind_reclassified = np.array(v_wind_list)
        u_wind_reclassified = np.array(u_wind_list)
        
        daily_data_dict[i].update({'v_wind':v_wind_reclassified, 'u_wind':u_wind_reclassified})
    
    fid.close()
    
    return daily_data_dict