# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:16:24 2020

@author: Jasper Dijkstra

Script contains function that are no longer used in TROPOMI plume detection algorithm:
    - GetTotalLandCells() - Get the total amount of grid cells above land
    - CheckCO() - Check CO concentration (in ppb) for given lat/lon combination
    - PlotDatabase() - Plot data from external databases (i.e. GFED, EDGAR, CAMS)

"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Local imports
import utilities as ut
import handling_EDGAR as edgar
import handling_GFED as gfed
import fetching_CAMS as cams
import masking_functions as mask




def GetTotalLandCells(boundaries, target_lat, target_lon):
    
    # Initialize array from data
    arr = np.ones((target_lat, target_lon))
    
    # Count amount of arrays that are land
    land = mask.land_sea_mask(arr, boundaries)
    total_land_cells = int(len(land[land != 0]))
    
    return total_land_cells


def checkCO(daily_data_dict, lat, lon):
    """
    Check CO concentration (in ppb) for given lat/lon combination

    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max and CO_ppb
    lat : float
        latitude.
    lon : float
        longitude.

    Returns
    -------
    CO_concentration : float
        Carbon Monoxide concentration (in ppb) as measured by TROPOMI on specified day

    """
    
    # Defining boundaries
    assert daily_data_dict['lat_min'] < lat < daily_data_dict['lat_max'] \
        and daily_data_dict['lon_min'] < lon < daily_data_dict['lon_max'], \
            'given lat/lon combination is out of reach!'
    
    lonrange = np.linspace(daily_data_dict['lon_min'], daily_data_dict['lon_max'], len(daily_data_dict['CO_ppb'][0]))
    latrange = np.linspace(daily_data_dict['lat_min'], daily_data_dict['lat_max'], len(daily_data_dict['CO_ppb']))
    
    # Get indices of lon and lat
    ilon = (np.abs(lonrange - lon)).argmin()
    ilat = (np.abs(latrange - lat)).argmin()
    CO_concentration = daily_data_dict['CO_ppb'][ilat, ilon]
    
    return CO_concentration



def PlotDatabase(dataset, in_path, out_path, bbox, xres, yres, day, month, year):
    """
    Function to plot data from the GFED, EDGAR and CAMS databases.
    
    Parameters
    ----------
    dataset : string
        Database from which data will be plotted. Can be 'GFED', 'EDGAR' or 'CAMS'
    in_path : string
        Path to folder containing file with data to be plotted.
    out_path : string
        Path to folder where output will be saved
    bbox : boundaries [lat_min, lat_max, lon_min, lon_max]
    day : day of TROPOMI overpass
    month : month of TROPOMI overpass
    year : year of TROPOMI overpass
    xres : resolution in lon direction
    yres : resolution in lat direction

    Returns
    -------
    Figure, stored in out_path
    
    """
    
    # Assert inputs are valid
    assert dataset in ['GFED', 'EDGAR', 'CAMS']
    if not os.path.isdir(in_path):
        print(f'Directory {in_path} was not found!')
        print(f'Creating {in_path} ...')
    in_path = ut.DefineAndCreateDirectory(in_path)
    out_path = ut.DefineAndCreateDirectory(out_path)
    
    # Check which database to read/open
    if dataset == 'GFED':
        data = gfed.OpenGFED(in_path, bbox, day, month, year, xres, yres)
        label = 'g C m^2 month^-1'
    if dataset == 'EDGAR':
        in_path = os.path.join(in_path + r'v432_CO_2012_IPCC_1A2.0.1x0.1.nc')
        data = edgar.OpenEDGAR(in_path, bbox, xres, yres)
        label = 'kg m-2 s-1'
    if dataset == 'CAMS':
        data = cams.FetchCams(in_path, dates=[day, month, year], bbox = bbox, xres = xres, yres = yres)
        label = 'ppb'
    
    # Mask data values equal to zero
    data = ma.array(data, mask=(data == 0))

    # Latitude and Longitudinal ranges
    latrange = np.linspace(bbox[0], bbox[1], data.shape[0])
    lonrange = np.linspace(bbox[2], bbox[3], data.shape[1])
    
    lon, lat = np.meshgrid(lonrange, latrange)
    
    # plot with cartopy
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # Add some cartopy features to the map
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
    ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3)
    
    cs = plt.pcolormesh(lon, lat, data, cmap='rainbow', transform=ccrs.PlateCarree())
    cbaxes = fig.add_axes([0.2, 0.03, 0.6, 0.03]) 
    cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
    cb.set_label(label)
    
    # Save the figure
    out_name = os.path.join(out_path + rf'{dataset}_{bbox[0]}{bbox[1]}{bbox[2]}{bbox[3]}.png')
    fig.savefig(out_name, bbox_inches='tight')#, dpi=1200)
    
    # Close the figure and clear the axes
    plt.cla()
    plt.clf()
    plt.close()
    
    return

