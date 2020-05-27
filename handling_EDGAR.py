# -*- coding: utf-8 -*-
"""
Created on Tue May 26 19:20:08 2020

@author: jaspd

Open EDGAR


Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). 
Emission Databasefor Global Atmospheric Research (EDGAR), "http://edgar.jrc.ec.europe.eu"

"""

from netCDF4 import Dataset
import numpy as np

# Local imports
import raster_tools as rast


def OpenEDGAR(path, bbox, xres, yres):
    # Open dataset
    fid = Dataset(path)
    
    # Retrieve latitude and longitudes
    lon = fid.variables['lon'][:]
    lat = fid.variables['lat'][:]
    CO_emissions = fid.variables['emi_co'][:] #kg m-2 s-1
    
    # Close dataset
    fid.close()
    
    # Get indeces of latitudes and logitudes where desired extent is True
    ilon = np.where((lon > bbox[2]) & (lon < bbox[3]))
    ilat = np.where((lat > bbox[0]) & (lat < bbox[1]))
    
    # Clip lon, lat and CO_emissions to desired extent 
    lon = lon[(lon > bbox[2]) & (lon < bbox[3])]
    lat = lat[(lat > bbox[0]) & (lat < bbox[1])]
    CO_emissions = CO_emissions[ilat[0], :]
    CO_emissions = CO_emissions[:, ilon[0]]
    
    # kg m-2 s-1 to ppb
    CO_emissions = CO_emissions / 366 # kg m-2 year (EDGAR year = 366 days)
    CO_emissions = CO_emissions * 1e+06 # ppb (1 kg per m3 = 1e+06 ppb)
    
    # Reclassify EDGAR data to desired extent
    CO_emissions = rast.ResampleArray(bbox, CO_emissions, xres, yres)
    
    #CO_emissions[CO_emissions > 0] = 1
    
    
    return CO_emissions
