# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:04:17 2020

@author: Jasper Dijkstra, edited from script from Todd Karin
"""

import h5py
import numpy as np
import os


def read_hdf5(day, month, year):

    global _lat
    global _lon
    global _mask
    global __file__
    
    # Load the data from file.
    #full_path = os.path.realpath('__file__')
    #path, _ = os.path.split(full_path)
    #
    #if not path.endswith(os.path.sep):
        #path += os.path.sep
    
    os.chdir(os.getcwd() + os.sep + 'GFED_fire_emissions_mask')
    updated_working_dir = os.getcwd()
    mask_filename = os.path.join(updated_working_dir, 'GFED4.1s_{}_beta.hdf5'.format(year))
    #mask_filename = os.path.join(path + r'/GFED_fire_emissions_mask/GFED4.1s_2017_beta.hdf5')
    
    #sys.path.append(path + r'/GFED_fire_emissions_mask')
    
    #mask_filename = 
    #mask_filename = os.path.join(r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\GFED_fire_emissions_mask\GFED4.1s_{}_beta.hdf5'.format(year))
    #mask_filename = os.path.join(full_path,r''.format(year))
    #mask_filename = os.path.join(path,'/GFED_fire_emissions_mask/GFED4.1s_{}_beta.hdf5'.format(year))
    
    f = h5py.File(mask_filename, mode='r')
        
    # Get daily (mask) and lat/lon data from hdf5 file
    _monthly_data = f['/emissions/{}/C'.format(month)][:]
    _daily_fraction = f['/emissions/{}/daily_fraction/day_{}'.format(month, day)][:]
    _mask = np.multiply(_monthly_data, _daily_fraction)
    _mask[_mask == 0] = np.NaN # Set all zero values to 'nan'
    _lat = f['/lat'][:]
    _lat = _lat[:, 0]
    _lon = f['/lon'][:]
    _lon = _lon[0]
    
    return


def is_wildfire(lat, lon, day, month, year):
    
    read_hdf5(day, month, year)
    
    # Convert latitude from TROPOMI to index on mask 
    lat = np.array(lat)
    if np.any(lat>90):
        raise ValueError('latitude must be <= 90')
    if np.any(lat<-90):
        raise ValueError('latitude must be >= -90')

    lat[lat > _lat.max()] = _lat.max()
    lat[lat < _lat.min()] = _lat.min()   
    
    #lat_i = ((lat - _lat[0])/(_lat[1]-_lat[0])).astype('int')
    lat_i = ((lat - _lat[0])/(_lat[1]-_lat[0])).astype('int')
    
    
    # Now convert longitude from TROPOMI to index on mask
    lon = np.array(lon)
    if np.any(lon > 180):
        raise ValueError('longitude must be <= 180')
    if np.any(lon < -180):
        raise ValueError('longitude must be >= -180')

    lon[lon > _lon.max()] = _lon.max()
    lon[lon < _lon.min()] = _lon.min()

    lon_i = ((lon - _lon[0]) / (_lon[1] - _lon[0])).astype('int')
    
    
    mask_array = np.logical_not(_mask[lat_i,lon_i])
    
    return mask_array

"""
#--------------------
# EXAMPLE
#--------------------



ary = is_wildfire(-50, lonrange, 13, 10, 2018)



read_hdf5(13, 10, 2018)

lat_min = -50
lat_max = 0
lon_min = 100
lon_max = 160
    
# Define timescope
day = 13
month = 10
year = 2018

# Defining the lat and longitudinal ranges
#latrange = np.linspace(lat_min, lat_max, len(daily_data[3]['CO_ppb']))
#lonrange = np.linspace(lon_min, lon_max, len(daily_data[3]['CO_ppb'][0]))


# Check for all lat/lon values if they correspond with wildfire emissions (1) or not (0)
mask = []
#array = daily_data[3]['CO_ppb']
for lat in range(len(latrange)):
    read_hdf5(day, month, year)
    is_wildfire = is_wildfire(latrange[lat], lonrange)
    is_wildfire = is_wildfire.astype(int)
    mask.append(is_wildfire)
mask = np.array(mask)

mask = np.multiply(mask, array)

"""
