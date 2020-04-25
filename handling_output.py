# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:42:35 2020

@author: Jasper Dijkstra edited from script from S. Houweling

This script contains functions to store the coordinates of plumes in a textfile,
and to create figures of plume locations ('mask') or atmospheric CO concentrations ('xCO')

"""

import numpy.ma as ma
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import utilities as ut


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
    
    # Get the indices of center of plume
    plume_mask = daily_data_dict['plume_mask']
    #indices = np.where(plume_mask == 1) 
    plumes = plume_mask > 0 # Define data to be labelled
    labels, nlabels = ndimage.label(plumes) # Label data
    indices = ndimage.center_of_mass(plume_mask, labels, np.arange(nlabels) + 1) # Define center of mass of labeled plumes
    max_xco = ndimage.measurements.maximum_position(field_t, labels, np.arange(nlabels) + 1)
    mean_xco = ndimage.measurements.maximum_position(field_t, labels, np.arange(nlabels) + 1)
    plume_size = ndimage.labeled_comprehension(plume_mask, labels, np.arange(nlabels) + 1, np.sum, float, 0)
    
    # Generate coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    # Write to txt file
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    curr_time = ut.GetCurrentTime()
    total = len(indices)
    filename = coord_directory + \
        'Plume_coordinates_{}_{}_{}.txt'.format(month, day, year)
    
    headerstring = """#----------------------------------------
#----------------------------------------
This file was automatically generated at: {}/{}/{} {}:{}

This file contains a list with information on Carbon Monoxide plumes at {}/{}/{}, between:
longitudes: [{}, {}] 
latitudes: [{}, {}] 

column descriptions:
- latitude:     Latitude of center of plume
- longitude:    Longitude of center of plume
- grid_cells:   Amount of grid cells (~7x7km) in plume
- CO_max:       Highest Carbon Monoxide concentration measured in plume (ppb)
- CO_average:   Average Carbon Monoxide concentration measured in plume (ppb)

Total amount of plumes identified: {}         
             
#----------------------------------------
#----------------------------------------
latitude, longitude, grid_cells, CO_max, CO_average,
""".format(curr_time['year'], curr_time['month'], curr_time['day'], \
    curr_time['hour'], curr_time['minute'], month, day, year, \
        lon_min, lon_max, lat_min, lat_max, total)
    
    f = open(filename, 'w+')
    f.write(headerstring)

    for i in range(len(indices)):
        x = int(indices[i][1])
        y = int(indices[i][0])
        x_co_max = int(max_xco[i][1])
        y_co_max = int(max_xco[i][0])
        f.write("{}, {}, {}, {}, {},\n".format(lat[y, x], lon[y, x], plume_size[i], field_t[y_co_max, x_co_max], mean_xco[i]))
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
        # Add coastlines
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
        ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3) 
        
        cs = plt.pcolormesh(lon, lat, field_mt, cmap='rainbow', transform=ccrs.PlateCarree())
        cbaxes = fig.add_axes([0.2, 0.1, 0.6, 0.03]) 
        cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
        cb.set_label('ppb')
        #ax.coastlines()
    
    elif figtype == 'mask':
        # Load topographical features from cartopy
        rivers_10m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '10m')
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
        ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m') 
        states_50m = cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces_lines','50m')
        lakes_50m = cfeature.NaturalEarthFeature('physical', 'lakes', '50m')
        
        # Add the topographical features to the map
        ax.add_feature(ocean_50m, edgecolor = 'face', facecolor = cfeature.COLORS['water'], zorder=1) 
        ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3)
        ax.add_feature(rivers_10m, facecolor='None',linewidth=0.25, edgecolor=cfeature.COLORS['water'],zorder=3)
        ax.add_feature(lakes_50m, edgecolor='k',linewidth=0.25,facecolor='None',zorder=3) 
        ax.add_feature(states_50m, edgecolor='gray',linewidth=0.25,facecolor='None',zorder=3)
        ax.add_feature(cfeature.BORDERS, edgecolor='#666666',linewidth=0.3,zorder=3)
        ax.patch.set_facecolor('None')
        
        # Creating the figure
        colors = ['white', 'red']
        cmap = ListedColormap(colors)
        cs = plt.pcolormesh(lon, lat, daily_data_dict['plume_mask'], cmap = cmap, transform=ccrs.PlateCarree())
        cbaxes = fig.add_axes([0.27, 0.05, 0.1, 0.03]) 
        cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal')
        cb.set_ticks([0, 1])
        cb.set_ticklabels(["no plume", "plume"])

    
    # Title in figure or not
    if title != None:
        plt.title(title)
    
    #ax.coastlines()
    plt.ioff() # Preventing figures from appearing as pop-up
    
    # Saving figure
    curr_time = ut.GetCurrentTime()
    fig.savefig(figue_directory + r'fig_{}_{}_{}_{}___{}{}{}{}{}{}.png'.format(figtype, month, day, year, \
            curr_time['year'], curr_time['month'], curr_time['day'], curr_time['hour'], curr_time['minute'], curr_time['second']), bbox_inches='tight')
    #print('figure saved at: {}'.format(saving_path))
    plt.close()
    
    return

