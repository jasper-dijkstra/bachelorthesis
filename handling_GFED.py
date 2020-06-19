# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:30:44 2020

@author: Jasper Dijkstra

Function to handle data from GFED on emissions (monthly, daily, 3-hourly) and burned area (with small fires).

Required data can be downloaded via:
    https://www.geo.vu.nl/~gwerf/GFED/GFED4/

"""

import h5py
import os
import numpy as np

# Local imports
import raster_tools as rast



def OpenGFED(path, bbox, day, month, year, xres, yres):
    """
    This function reads and resamples HDF5 data from the Global Fire Emissions Database (GFED),
    GFED data from 2017 onwards, is only available as Beta version. 
    Therefore this function makes no distinction in plume intensity.
    
    Parameters
    ----------
    path : string
        folder where gfed data is stored.
    bbox : list, tuple
        list, tuple with desired data extents [lat_min, lat_max, lon_min, lon_max]
    day, month, year : int, str
        desired day, month and year GFED data needs to be acquired from
    xres, yres : int, float
        desired resolution in longitudal direction (xres) and latitudal direction (y)
        
    Returns
    -------
    2d np.ndarray of desired resolution with GFED data for specified time scope (in g C m-2 month-1)
    
    """

    # Read the data from the hdf5 file
    filename = str(os.path.join(path + rf'GFED4.1s_{year}_beta.hdf5'))
    f = h5py.File(filename, mode='r')
    
    # Get daily fire and lat/lon data from hdf5 file
    monthly_data = f['/emissions/{}/C'.format(month)][:]
    daily_fraction = f['/emissions/{}/daily_fraction/day_{}'.format(month, day)][:]
    data = np.multiply(monthly_data, daily_fraction)
    lat = f['/lat'][:,0]
    lon = f['/lon'][0,:]
    
    # Close the hdf5 file to save memory
    f.close()
    
    # Get indeces of latitudes and logitudes where desired extent is True
    ilon = np.where((lon > bbox[2]) & (lon < bbox[3]))
    ilat = np.where((lat > bbox[0]) & (lat < bbox[1]))
    
    # Clip lon, lat and CO_emissions to desired extent 
    lon = lon[(lon > bbox[2]) & (lon < bbox[3])]
    lat = lat[(lat > bbox[0]) & (lat < bbox[1])]
    data = data[ilat[0], :]
    data = data[:, ilon[0]]
    
    # Resample data to desired resolution
    CO_emissions = rast.ResampleArray(bbox, data, xres, yres)
    
    # Flipping the data, as for some reason it is upside down
    CO_emissions = np.flipud(CO_emissions)
    
    return CO_emissions