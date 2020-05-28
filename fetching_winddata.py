#!/usr/bin/env python
"""
Created on Mon Apr 27 21:38:48 2020

@author: Jasper Dijkstra

This script contains functions to:
    1. Download ERA5 meteodata on u- and v-winds for a desired atmospheric pressure
        level, time- and aerial scope.
    2. Process the data obtained with function 1, to np arrays the size of input
    
"""

from datetime import datetime, timedelta
from netCDF4 import Dataset
import cdsapi
import numpy as np
from scipy.ndimage import gaussian_filter
import os

# Local imports
import utilities as ut
import raster_tools as raster


def GetTimeZoneExtent(lon_min, lon_max, lat_min, lat_max, UTCzones, hour):
    # Get the extents of the timezone in which data has to be requested
    
    # Create Meshgrid with current coordinates
    lonrange = np.linspace(lon_min, lon_max, len(UTCzones[0]))
    latrange = np.linspace(lat_min, lat_max, len(UTCzones))
    lon, lat = np.meshgrid(lonrange, latrange)
    
    # Get longitude (x) coorinates
    lon[UTCzones != hour] = np.nan # Set all lon's that are not <hour>, to 'nan' 
    xmin = np.nanmin(lon)
    xmax = np.nanmax(lon)
    
    # Get latitude (y) coordinates
    #lat[UTCzones != hour] = np.nan  # Set all lat's that are not <hour>, to 'nan' 
    ymin = np.nanmin(lat)
    ymax = np.nanmax(lat)
    

    # Make sure extents fall within boundary
    if not -180 <= xmin <= 180 or not -180 <= xmax <= 180 or not xmin <= xmax:
        print('maximum longitude cannot be smaller than or equal to minimum, and should be within range (-180, 180)! Not able to place request!')
        return 'Error'
    if not -90 <= ymin <= 90 or not -90 <= ymax <= 90 or not ymin <= ymax:
        print('maximum latitude cannot be smaller than or equal to minimum, and should be within range (-90, 90)! Not able to place request!')
        return 'Error'
    
    
    zone_extents = {'hour':hour, \
                    'lon_min':xmin, 'lon_max':xmax, \
                        'lat_min':ymin, 'lat_max':ymax}
    
    return zone_extents  

    
def GetCorrectDates(zone_extents, daily_data_dict, timerange):
    day = str(daily_data_dict['day'])
    month = str(daily_data_dict['month'])
    year = str(daily_data_dict['year'])
    hour = zone_extents['hour']
    
    capture_time = datetime.strptime(f'{year} {month} {day} {hour}', '%Y %m %d %H')
    meteodata_range = capture_time - timedelta(hours = timerange-1)
    
    if zone_extents['hour'] - timerange + 1 < 0:
        hours_today = np.linspace(0, zone_extents['hour'], zone_extents['hour']+1)
        hours_left = abs(zone_extents['hour'] - timerange + 1)
        hours_yesterday = np.linspace(24-hours_left, 23, hours_left)
        hours_raw = list(hours_yesterday) + list(hours_today)
    else:
        hours_raw = list(np.linspace(zone_extents['hour']-timerange+1, zone_extents['hour'], timerange))
    
    hours = list()
    for h in hours_raw:
        if h < 10:
            hours.append(str(int(h)).zfill(2) + ':00')
        else:
            hours.append(str(int(h)) + ':00')
    
    days = [meteodata_range.day, capture_time.day]
    months = [meteodata_range.month, capture_time.month]
    years = [meteodata_range.year, capture_time.year]
    
    # Remove duplicates in list
    for list_ in [days, months, years]:
        if list_[0] == list_[1]:
            list_.pop(0)

    return [years, months, days, hours]


def DownloadMeteo(zone_extents, timescope, basepath, pressure_level=850):
    """
    This function downloads ERA5 data making use of the CDS API from Copernicus,
    as netCDF4 files for the time- and aerial scope defined in daily_data_dict
    
    Parameters
    ----------
    zone_extents : output of GetTimeZoneExtent()
    timescope : output of GetCorrectDates()
    pressure_level : int
        Atmospheric pressure level (hPa) at which data has to be downloaded
    basepath : string
        Path to folder where meteodata will be stored
        
        
    Returns
    -------
    downloaded .nc file in: (<working_dir>/store_meteo_data/<subdir>).

    """
     
    # redefining earlier results
    lat_min = zone_extents['lat_min']
    lat_max = zone_extents['lat_max']
    lon_min = zone_extents['lon_min']
    lon_max = zone_extents['lon_max']
    
    years = timescope[0]
    months = timescope[1]
    days = timescope[2]
    hours = timescope[3]
    
    # Setting output directory
    #curr_directory = os.getcwd()
    out_dir = os.path.join(basepath, rf'03_Store_meteo_data/{years[0]}/{months[0]}/{days[0]}')
    ut.DefineAndCreateDirectory(out_dir)
    filename = f'ERA5_{hours[-1][:2]}h_Lon[{lon_min}_{lon_max}]_Lat[{lat_min}_{lat_max}]_{pressure_level}hPa.nc'
    
    if os.path.isfile(os.path.join(out_dir, filename)) == True:
        print('ERA5 meteodata already exists!')
        print('No new attempt to download the data will be undertaken!')
        return os.path.join(out_dir, filename)
    
    # Actual request:
    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'area': [lat_min, lon_min, lat_max, lon_max],
            'time': hours,
            'day': days,
            'month': months,
            'year': years,
            'pressure_level': str(pressure_level),
            'variable': ['u_component_of_wind', 'v_component_of_wind'],
        },
        os.path.join(out_dir, filename))
    
    return os.path.join(out_dir, filename)


def OpenERA5(path, zone_extents, timerange):
    """
    Function to read meteodata opened with DownloadMeteo().
    Returns: [lon_meteo, lat_meteo, u_wind, v_wind]

    """
    # Open dataset
    fid = Dataset(path)
    
    lon_meteo = fid.variables['longitude'][:]
    lat_meteo = fid.variables['latitude'][:]
    time = fid.variables['time'][:]
    
    # Check if data is spread over multiple days 
    # (if so, make sure only subsequent hours are taken into account) 
    if len(time) != timerange:
        # Define the upper and lower indices of subsequent timerange
        lower_indices = zone_extents['hour'] + 1
        upper_indices = lower_indices + timerange 
        
        # Take the mean over all hours for u and v wind
        u_wind = fid.variables['u'][lower_indices:upper_indices]
        u_wind = np.mean(u_wind, axis=0)
        v_wind = fid.variables['v'][lower_indices:upper_indices]
        v_wind = np.mean(v_wind, axis=0)
        
    else:
        # Take the mean over all hours for u and v wind
        u_wind = fid.variables['u'][:]
        u_wind = np.mean(u_wind, axis=0)
        v_wind = fid.variables['v'][:]
        v_wind = np.mean(v_wind, axis=0)
    
    fid.close()
    
    return [lon_meteo, lat_meteo, u_wind, v_wind]


def FindNearest(array, value):
    """
    Find index of <array> value that lies closest to input <value> 
    """
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    return idx


def RescaleWindArray(daily_data_dict, zone_extents, ERA5):
    """
    Rescale wind data to resolution of TROPOMI data

    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>].
    ERA5 : output of OpenERA5()

    Returns
    -------
    uwind : np.array with requested ERA5 data on uwind on TROPOMI resolution
    vwind : np.array with requested ERA5 data on vwind on TROPOMI resolution
    lon_m_range : longitudinal range of ERA5 data on wind 
    lat_m_range : latitudal range of ERA5 data on wind 

    """
    # Redefine ERA5 data
    lon_meteo = ERA5[0]
    lat_meteo = ERA5[1]
    u_wind = ERA5[2]
    v_wind = ERA5[3]
    
    # Defining target boundaries from zone_extents
    lat_min = zone_extents['lat_min']
    lat_max = zone_extents['lat_max']
    lon_min = zone_extents['lon_min']
    lon_max = zone_extents['lon_max']
    bbox = [lat_min, lat_max, lon_min, lon_max]
    
    # Calculate array resolution (km) for lat and lon direction
    lon_res_out = 110*((daily_data_dict['lon_max']-daily_data_dict['lon_min'])/daily_data_dict['target_lon'])
    lat_res_out = 110*((daily_data_dict['lat_max']-daily_data_dict['lat_min'])/daily_data_dict['target_lat'])
    
    uwind = raster.ResampleArray(bbox, u_wind, lon_res_out, lat_res_out)
    vwind = raster.ResampleArray(bbox, v_wind, lon_res_out, lat_res_out)

    lon_m_range = np.linspace(np.min(lon_meteo), np.max(lon_meteo), uwind.shape[1])
    lat_m_range = np.linspace(np.min(lat_meteo), np.max(lat_meteo), uwind.shape[0])
    lon_m, lat_m = np.meshgrid(lon_m_range, lat_m_range)
    
    # Flatten the arrays to 1D
    uwind = uwind.flatten()
    vwind = vwind.flatten()
    lon_m = lon_m.flatten()
    lat_m = lat_m.flatten()
    
    return uwind, vwind, lon_m, lat_m


def FetchWindData(daily_data_dict, pressure, timerange, basepath):
    """
    Function to fetch meteorological data from ECMWF ERA5 on uwind and vwind:
        - for specified air pressure (in hPa)
        - for specified amount of hours (<timerange>):
            hours that will be taken into account before capture of TROPOMI data

    Parameters
    ----------
    daily_data_dict : TYPE
        DESCRIPTION.
    pressure : integer
        Atmospheric pressure in hPa.
    timerange : int
        Amount of hours to take into account before TROPOMI observation was done.
    basepath : string
        Path to folder where meteodata will be stored

    Returns
    -------
    daily_data_dict : dictionary
        Same as input, only u_wind and v_wind have been added

    """
    
    # Redefine output extents
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Initialize arrays for lat/lon
    lonrange = np.linspace(lon_min, lon_max, daily_data_dict['target_lon'])
    latrange = np.linspace(lat_min, lat_max, daily_data_dict['target_lat'])
    lon, lat = np.meshgrid(lonrange, latrange)
    
    
    # Redefine array with hours of TROPOMI data catpure to UTCzones
    UTCzones = daily_data_dict['hour']
    
    # Create a list with all hours within UTCzones
    flat = UTCzones.flatten() # Change the UTCzones 2D array to one dimension
    flat = flat[~np.isnan(flat)].astype(int) # Remove all 'nan' values
    count = np.bincount(flat) # Count the amount of values in hour array
    treshold = np.quantile(count, 0.5) # Limit nuisances by removing not frequent hours
    
    hourslist = list() # Initialize hourslist
    for i in range(len(count)):
        if count[i] > treshold: # Make sure value frequency is above treshold
            hourslist.append(i)
    
    # Initialize arrays for output wind arrays
    u_wind = np.zeros(UTCzones.shape)
    v_wind = np.zeros(UTCzones.shape)
    count_arr = np.zeros(UTCzones.shape)
    
    
    for h in range(0, len(hourslist)):
        # 1. Get extents of TROPOMI orbit
        zone_extents = GetTimeZoneExtent(lon_min, lon_max, lat_min, lat_max, UTCzones, hour=hourslist[h])
        if zone_extents == 'Error':
            continue
        
        # 2. Get timerange in which meteodata will be requested
        timescope = GetCorrectDates(zone_extents, daily_data_dict, timerange=int(timerange))
        
        # 3. Request meteodata via API
        path = DownloadMeteo(zone_extents, timescope, basepath, pressure_level=pressure)
        
        # 4. Read requested meteodata as np.arrays
        ERA5 = OpenERA5(path, zone_extents, timerange=int(timerange))
        
        # 5. Rescale meteodata to TROPOMI resolution
        u_wind_m, v_wind_m, lon_m, lat_m = RescaleWindArray(daily_data_dict, zone_extents, ERA5)
        
        # 6. Fit meteodata in TROPOMI array
        
        u = np.zeros(UTCzones.shape)
        v = np.zeros(UTCzones.shape)
        c = np.zeros(UTCzones.shape)
        
        for iobs in range(len(lon_m)):
            
            # Find the index in output that is closest to input
            ilon = FindNearest(lonrange, lon_m[iobs])
            ilat = FindNearest(latrange, lat_m[iobs])
            
            # Append to hourly u_wind, v_wind & count frame
            u[ilat, ilon] += u_wind_m[iobs]
            v[ilat, ilon] += v_wind_m[iobs]
            c[ilat, ilon] += 1
            
        mask = (UTCzones != h) | (UTCzones == np.nan)
        u_wind += (u * mask)
        v_wind += (v * mask)
        count_arr += (c * mask)
        
    idx = (count_arr > 0)
    u_wind[idx] = u_wind[idx]/count_arr[idx]
    v_wind[idx] = v_wind[idx]/count_arr[idx]
    
    # Apply smoothing filter here?
    #u_wind = gaussian_filter(u_wind, 6)
    #v_wind = gaussian_filter(v_wind, 6)
    
    return u_wind, v_wind

