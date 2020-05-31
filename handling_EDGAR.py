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
import numpy.ma as ma
from scipy import constants

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
    
    # Only highest emissions are relevant
    nanCO = np.copy(CO_emissions)
    nanCO[nanCO == 0] = np.nan
    mean = np.nanmean(nanCO)
    stdev = np.nanstd(nanCO)
    CO_emissions[CO_emissions < mean+stdev] = 0
    
    # Get indeces of latitudes and logitudes where desired extent is True
    ilon = np.where((lon > bbox[2]) & (lon < bbox[3]))
    ilat = np.where((lat > bbox[0]) & (lat < bbox[1]))
    
    # Clip lon, lat and CO_emissions to desired extent 
    lon = lon[(lon > bbox[2]) & (lon < bbox[3])]
    lat = lat[(lat > bbox[0]) & (lat < bbox[1])]
    CO_emissions = CO_emissions[ilat[0], :]
    CO_emissions = CO_emissions[:, ilon[0]]
    
# =============================================================================
#     # kg m-2 s-1 to kg m-2 day-1
#     average = 85 / 1000000 #ppb multiplied to kg/m3
#     CO_emissions = CO_emissions * 3600 * 24 # kg m-2 day-1 (EDGAR year = 366 days)
#     CO_emissions = CO_emissions * (average/11000)
#     CO_emissions = CO_emissions * 1e+06
# =============================================================================

    # Reclassify EDGAR data to desired extent
    CO_emissions = rast.ResampleArray(bbox, CO_emissions, xres, yres)
    
    # As we are only interested to know if emissions happened, return emissions true(1)/false(0)
    CO_emissions[CO_emissions > 0] = 1
    
    
    return CO_emissions


