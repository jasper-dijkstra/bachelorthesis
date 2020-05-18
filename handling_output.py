# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 11:42:35 2020

@author: Jasper Dijkstra edited from script from S. Houweling

This script contains functions to store the coordinates of plumes in a textfile,
and to create discrete or continuous georeferenced figures of np.arrays

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
    mean_xco = ndimage.measurements.mean(field_t, labels, np.arange(nlabels) + 1)
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
    
    headerstring = f"""#----------------------------------------
#----------------------------------------
This file was automatically generated at: {curr_time['year']}/{curr_time['month']}/{curr_time['day']} {curr_time['hour']}:{curr_time['minute']}

This file contains a list with information on Carbon Monoxide plumes at {month}/{day}/{year}, between:
longitudes: [{lon_min}, {lon_max}] 
latitudes: [{lat_min}, {lat_max}] 

column descriptions:
- latitude:     Latitude of center of plume
- longitude:    Longitude of center of plume
- grid_cells:   Amount of grid cells (~7x7km) in plume
- CO_max:       Highest Carbon Monoxide concentration measured in plume (ppb)
- CO_average:   Average Carbon Monoxide concentration measured in plume (ppb)

Total amount of plumes identified: {total}         
             
#----------------------------------------
#----------------------------------------
latitude, longitude, grid_cells, CO_max, CO_average,
"""
    
    f = open(filename, 'w+')
    f.write(headerstring)

    for i in range(len(indices)):
        x = int(indices[i][1])
        y = int(indices[i][0])
        x_co_max = int(max_xco[i][1])
        y_co_max = int(max_xco[i][0])
        f.write("{}; {}; {}; {}; {};\n".format(lat[y, x], lon[y, x], plume_size[i], field_t[y_co_max, x_co_max], mean_xco[i]))
    f.close()
    
    return



def CreateTargetLatLonGrid(daily_data_dict, figtype):
    # Defining boundaries
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Setting target lon- and latitude
    nlon_t = len(daily_data_dict[figtype][0])
    nlat_t = len(daily_data_dict[figtype])
    
    # Generate coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    return lon, lat


def Masking(masking, daily_data_dict, figtype):
    if masking:
        count_t = daily_data_dict['count_t']
        mask = (count_t == 0) # Masking all zero values in count_t array
        field_mt = ma.array(daily_data_dict[figtype], mask=mask) # Apply this mask to figtype as well
    else:
        field_mt = daily_data_dict[figtype]
    
    return field_mt


def Title(ax, title, figtype, year, month, day):
    if title != None:
        ax.set_title(title, horizontalalignment='center', verticalalignment='top')
    else:
        ax.set_title(f'Map of {figtype}, at: {year}/{month}/{day}', horizontalalignment='center', verticalalignment='top')
    
    return ax


def ExportFig(fig, figure_directory, figtype, month, day, year):
    curr_time = ut.GetCurrentTime()
    fig.savefig(figure_directory + rf"fig_{figtype}_{month}_{day}_{year}___{curr_time['year']}{curr_time['month']}{curr_time['day']}{curr_time['hour']}{curr_time['minute']}{curr_time['second']}.png", bbox_inches='tight')
    
    return

def CreateMaskMap(daily_data_dict, figtype, figure_directory, title=None, \
                  labels=["no plume", "plume"], colors=['white', 'red']):
    """
    
    Fucntion to create a discrete color map of a 2D np.array
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max, day, month, year,
        and np.array with values to plot and count_t.
    figtype : string
        Name as in daily_data_dict of the np.array containing values to plot.
    figure_directory : string
        Directory where figures will be stored.
    labels : list with strings, optional
        Labels for each tick on colorbar. The default is ["no plume", "plume"].
        NOTE! 'labels' must have the same length as 'colors'!
    colors : list with strings, optional
        colors for each tick on colorbar. The default is ['white', 'red'].
        NOTE! 'colors' must have the same length as 'labels'!

    Returns
    -------
    Map saved as png in figure_directory.

    """
    # Check if labels and colors have got the same lengths
    if not len(labels) == len(colors):
        print("Parameters 'labels' and 'colors' must have the same length!")
        print("Setting to default: labels=[no plume, plume], colors=[white, red]")
    
    # Retrieving month, day and year
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    
    # Create a longitude and latitude np.meshgrid in target resolution
    lon, lat = CreateTargetLatLonGrid(daily_data_dict, figtype)
    
    # Create cartopy plot
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
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
    cmap = ListedColormap(colors)
    cs = plt.pcolormesh(lon, lat, daily_data_dict[figtype], cmap = cmap, transform=ccrs.PlateCarree())
    cbaxes = fig.add_axes([0.27, 0.05, 0.1, 0.03]) 
    cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal')
    cb.set_ticks([0, len(colors)-1])
    cb.set_ticklabels(labels)

    # Set title of figure
    ax = Title(ax, title, figtype, year, month, day)
    
    plt.ioff() # Preventing figures from appearing as pop-up
    
    # Saving figure
    ExportFig(fig, figure_directory, figtype, month, day, year)

    plt.close()


def CreateColorMap(daily_data_dict, figtype, figure_directory, masking=True, \
                   labeltag='values', title=None):
    """
    
    Fucntion to create a continuous color map of a 2D np.array
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max, day, month, year,
        and np.array with values to plot and count_t.
    figtype : string
        Name as in daily_data_dict of the np.array containing values to plot.
    figure_directory : string
        Directory where figures will be stored.
    masking : bool, optional
        DESCRIPTION. The default is True.
    labeltag : string, optional
        Tag for the colorbar of the plot (for example 'ppb' or '(g C / m^2)'.
        The default is 'values'.
    title : string, optional
        Title to be given to the plot. The default is None.

    Returns
    -------
    Map saved as png in figure_directory.

    """

    
    # Retrieving month, day and year
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    
    # Create a longitude and latitude np.meshgrid in target resolution
    lon, lat = CreateTargetLatLonGrid(daily_data_dict, figtype)
    
    # plot with cartopy
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Apply mask
    field_mt = Masking(masking, daily_data_dict, figtype)
    
    # Add some cartopy features to the map
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
    ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3) 
        
    cs = plt.pcolormesh(lon, lat, field_mt, cmap='rainbow', transform=ccrs.PlateCarree())
    cbaxes = fig.add_axes([0.2, 0.1, 0.6, 0.03]) 
    cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
    cb.set_label(labeltag)
    
    # Set title of figure
    ax = Title(ax, title, figtype, year, month, day)
    
    plt.ioff() # Preventing figures from appearing as pop-up
    
    # Saving figure
    curr_time = ut.GetCurrentTime()
    fig.savefig(figure_directory + r'fig_{}_{}_{}_{}___{}{}{}{}{}{}.png'.format(figtype, month, day, year, \
            curr_time['year'], curr_time['month'], curr_time['day'], curr_time['hour'], curr_time['minute'], curr_time['second']), bbox_inches='tight')
    plt.close()  



def CreateWindVector(daily_data_dict, figtype, figure_directory, masking=False, \
                   labeltag='m/s', title=None, skip=30, vector_only=False):
        
    """
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max, day, month, year,
        u_wind and v_wind.
    figure_directory : string
        Directory where figures will be stored.
    labeltag : string, optional
        Tag for the colorbar of the plot. The default is 'm/s'.
    title : string, optional
        Title to be given to the plot. The default is None.
    vector_only : bool, optional
        Only draw wind vectors if True. If False also colored wind speed is drawn.
        The default is False.
    skip : int, optional
        Every <skip> arrow is drawn in figure. The default is 30.

    Returns
    -------
    Map saved as png in figure_directory.

    """
    
    # Retrieving month, day and year
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    
    # Define U and V wind 
    U = daily_data_dict['u_wind']
    V = daily_data_dict['v_wind']
    
    windspeed = (U ** 2 + V ** 2) ** 0.5
    #direction = 180 + (180/np.pi)*np.arctan2(V,U)
    
    # Create a longitude and latitude np.meshgrid in target resolution
    lon, lat = CreateTargetLatLonGrid(daily_data_dict, 'u_wind')
    
    # plot with cartopy
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Add some cartopy features to the map
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
    ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3) 

    # Draw vectors
    cs = plt.quiver(lon[::skip, ::skip], lat[::skip, ::skip], U[::skip, ::skip], V[::skip, ::skip], transform=ccrs.PlateCarree(), pivot = 'mid', zorder=2)
    
    # Color wind speed as background as well
    if not vector_only:
        if figtype == 'windspeed':
            cs = plt.pcolormesh(lon, lat, windspeed, cmap='rainbow', transform=ccrs.PlateCarree(), zorder = 1)
        else:
            # Apply mask
            field_mt = Masking(masking, daily_data_dict, figtype)
            cs = plt.pcolormesh(lon, lat, field_mt, cmap='rainbow', transform=ccrs.PlateCarree(), zorder = 1)
        cbaxes = fig.add_axes([0.2, 0.1, 0.6, 0.03]) 
        cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
        cb.set_label(labeltag)
    
    # Set title of figure
    ax = Title(ax, title, figtype, year, month, day)
    
    plt.ioff() # Preventing figures from appearing as pop-up
    
    # Saving figure
    ExportFig(fig, figure_directory, figtype, month, day, year)
    
    plt.close() 

    return


def scatterplot(x, y, figure_directory, title=None, popup=False):
    plt.xlabel('buffersize')
    plt.ylabel('amount of plume grid cells')
    
    if title != None:
        plt.title(title)
    
    plt.scatter(x, y)
    curr_time = ut.GetCurrentTime()
    plt.savefig(figure_directory + r'scatter_{}{}{}{}{}{}.png'.format(curr_time['year'], curr_time['month'], curr_time['day'], curr_time['hour'], curr_time['minute'], curr_time['second']), bbox_inches='tight')
    plt.close()
    
    return