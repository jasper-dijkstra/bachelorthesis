# -*- coding: utf-8 -*-
"""
Created on Tue May 26 19:20:08 2020

@author: Jasper Dijkstra

Open the EDGAR database and only keep highest emissions as array

Source: European Commission, Joint Research Centre (JRC)/Netherlands Environmental Assessment Agency (PBL). 
Emission Databasefor Global Atmospheric Research (EDGAR), "http://edgar.jrc.ec.europe.eu"

"""

from netCDF4 import Dataset
import numpy as np

# Local imports
import raster_tools as rast


def OpenEDGAR(path, bbox, xres, yres):
    """
    This function reads and resamples data from the Emission Database for Global Atmospheric Research (EDGAR)
    
    Only highest emissions will be considered as having a significant influence on the column mixing ratio as observed by TROPOMI
    
    Parameters
    ----------
    path : string
        folder where gfed data is stored.
    bbox : list, tuple
        list, tuple with desired data extents [lat_min, lat_max, lon_min, lon_max]
    xres, yres : int, float
        desired resolution in longitudal direction (xres) and latitudal direction (y)
        
    Returns
    -------
    2d np.ndarray of desired resolution with EDGAR data for specified time scope (in kg m-2 s-1)
    
    """
    # EDGAR longitudes range from 0-360, while TROPOMI does -180-180
    # Therefore apply correction to bbox    
    edgar_bbox = [bbox[0], bbox[1], bbox[2] + 180, bbox[3] + 180]

    # Open dataset
    fid = Dataset(path)
    
    # Retrieve latitude and longitudes
    lon = fid.variables['lon'][:]
    lat = fid.variables['lat'][:]
    CO_emissions = fid.variables['emi_co'][:] #kg m-2 s-1
    
    # Put longitudes in correct order: [180 - 360] block to [-180 - 0]
    west = np.squeeze(CO_emissions[:,np.where(lon > 180)])
    east = np.squeeze(CO_emissions[:,np.where(lon <= 180)])
    CO_emissions = np.concatenate((west,east),axis=1)
    
    # Close dataset
    fid.close()
    
    # Only highest emissions are relevant to keep
    nanCO = np.copy(CO_emissions)
    nanCO[nanCO == 0] = np.nan
    mean = np.nanmean(nanCO)
    stdev = np.nanstd(nanCO)
    CO_emissions[CO_emissions < mean+stdev] = 0
    
    # Get indices of latitudes and logitudes where desired extent is True
    ilon = np.where((lon > edgar_bbox[2]) & (lon < edgar_bbox[3]))
    ilat = np.where((lat > edgar_bbox[0]) & (lat < edgar_bbox[1]))
    
    # Clip lon, lat and CO_emissions to desired extent 
    lon = lon[ilon[0]]
    lat = lat[ilat[0]]
    CO_emissions = CO_emissions[ilat[0], :]
    CO_emissions = CO_emissions[:, ilon[0]]

    # Reclassify EDGAR data to desired extent
    CO_emissions = rast.ResampleArray(bbox, CO_emissions, xres, yres)
    
    return CO_emissions
