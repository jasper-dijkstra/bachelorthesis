# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:42:35 2020

@author: Jasper Dijkstra edited from script from S. Houweling

This script contains functions to store the coordinates of plumes in a textfile,
and to create figures of plume locations ('mask') or atmospheric CO concentrations ('xCO')

"""

import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs



def NotePlumeCoordinates(daily_data_dict, coord_directory):
    """
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
    coord_directory : string
        directory where list (txt file) of plume coordinates will be stored.

    Returns
    -------
    Saves .txt file in coord_directory.

    """
    
    # Defining boundaries
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Deciding on the nlon_t and nlat_t
    field_t = daily_data_dict['CO_ppb']
    nlon_t = len(field_t[0])
    nlat_t = len(field_t)
    
    # Check the indices where a plume is detected
    plume_mask = daily_data_dict['plume_mask']
    indices = np.where(plume_mask == 1)
    #mask = (plume_mask == 0)
    #field_mt = ma.array(field_t, mask=mask)
    
    # Generate coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    # Write to txt file
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    total = len(indices[0])
    filename = coord_directory + 'Plume_coordinates_{}_{}_{}.txt'.format(month, day, year)
    
    headerstring = """#----------------------------------------
#----------------------------------------
This file contains a list with coordinates of plumes on {}/{}/{}, between:
longitudes: [{}, {}] 
latitudes: [{}, {}] 

All listed coordinates give the center of a gridcel of ~7x7km

Total amount of gridcells identified as plume: {}          
             
#----------------------------------------
#----------------------------------------
""".format(month, day, year, lon_min, lon_max, lat_min, lat_max, total)
    
    f = open(filename, 'w+')
    f.write(headerstring)

    for i in range(len(indices[0])):
        x = indices[1][i]
        y = indices[0][i]
        f.write("\nLat: {}, Lon: {}, xCO: {} ppb".format(lat[y, x], lon[y, x], field_t[y, x]))
    f.close()
    
    return




def CreateFigue(daily_data_dict, figue_directory, figtype, title=None):
    """
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
    figure_directory : string
        directory where figures will be stored.
    figtype : string
        'mask' or 'xCO' depending on the type of figure you want to make.
    title : string, optional
        title of the figure. The default is None.

    Returns
    -------
    saves figure in figure_directory.

    """
    
    if figtype not in ['mask', 'xCO']:
        print('figtype is not recognised, figure could not be created')
        # Add logging message
        return
    
    # Defining boundaries
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Retrieving month, day and year
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    
    # Deciding on the nlon_t and nlat_t
    field_t = daily_data_dict['CO_ppb']
    nlon_t = len(field_t[0])
    nlat_t = len(field_t)
    
    # make new masked array of valid averages
    count_t = daily_data_dict['count_t']
    mask = (count_t == 0)
    field_mt = ma.array(field_t, mask=mask)
    
    # Generate coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    # plot with cartopy
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Deciding what data will be plotted
    if figtype == 'xCO':
        cs = plt.pcolormesh(lon, lat, field_mt, cmap='rainbow', transform=ccrs.PlateCarree())
        cbaxes = fig.add_axes([0.2, 0.1, 0.6, 0.03]) 
        cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
        cb.set_label('ppb')
    
    # IMPROVE VISUALISATION OF MASK!
    elif figtype == 'mask':
        cs = plt.pcolormesh(lon, lat, daily_data_dict['plume_mask'], cmap = plt.cm.get_cmap('Reds', 2), transform=ccrs.PlateCarree())
        cbaxes = fig.add_axes([0.8, 0.75, 0.03, 0.05]) 
        cb = plt.colorbar(cs, cax = cbaxes, orientation = 'vertical', label = 'label')
        cb.set_label('no plume (0), plume (1)')
    
    # Title in figure or not
    if title != None:
        plt.title(title)
    
    ax.coastlines()
    plt.ioff() # Preventing figures from appearing as pop-up
    
    # Saving figure
    fig.savefig(figue_directory + r'fig_{}_{}_{}_{}.png'.format(figtype, month, day, year), bbox_inches='tight')
    #print('figure saved at: {}'.format(saving_path))
    plt.close()
    
    return

