# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 15:17:33 2020

@author: Jasper Dijkstra, based on a script from Mats Riet

Script to generate figures as subplots or seperate files for:
    - Total column CO as observed by TROPOMI
    - Plumes identified by the algorithm compared to the GFED and EDGAR databases
    - Wind data plotted as vector field

"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Local imports
import utilities as ut


def PlotFigures(daily_data_dict, basepath, subplots=True, mask_type = 'plume_mask'):
    """
    Plot figures of:
        - CO concentration as detected by TROPOMI
        - Plumes detected in the TROPOMI data, compared to GFED and EDGAR databases


    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max, day, month, year,
        and np.array with values to plot and count_t.
    basepath : str
        String to folder where output files will be stored.
    subplots : BOOL, optional
        If True, both maps are plotted in one figure (subplots)
        If False, two separate outputs are generated
        The default is True.

    Returns
    -------
    Figure(s) stored as png in '...'basepath'.../04_output/plume_figures'
    

    """
    
    # Assert input is valid
    assert mask_type in ['plume_mask', 'plumes_explained', 'plumes_incl_gfed_edgar'], 'Could not identify the plumes to plot'
    
    # 1: DEFINING OUTPUT NAME AND STORAGE LOCATION
    # Creating output directory
    out_dir = ut.DefineAndCreateDirectory(os.path.join(basepath + os.path.sep + r'04_output\plume_figures'))
    
    # Creating output filename
    curr_time = ut.GetCurrentTime()
    tropomi_datestamp = f"{daily_data_dict['month']}_{daily_data_dict['day']}_{daily_data_dict['year']}"
    current_timestamp = f"{curr_time['year']}{curr_time['month']}{curr_time['day']}{curr_time['hour']}{curr_time['minute']}{curr_time['second']}"
    out_name_mask = out_dir + rf"fig_plumes_{tropomi_datestamp}___{current_timestamp}.png"
    out_name_co = out_dir + rf"fig_CO_ppb_{tropomi_datestamp}___{current_timestamp}.png"
    out_name = out_dir + rf"fig_subplots_{tropomi_datestamp}___{current_timestamp}.png"

    
    # 2: CREATE LAT LON GRID TO GEOREFERENCE ALL GRID CELLS
    # Setting target lon- and latitude
    nlat_t, nlon_t = daily_data_dict['CO_ppb'].shape

    # Creating lat/lon meshgrid, to correctly place grid cells
    lon_t = np.linspace(daily_data_dict['lon_min'], daily_data_dict['lon_max'], nlon_t)
    lat_t = np.linspace(daily_data_dict['lat_min'], daily_data_dict['lat_max'], nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    # 3: PREPARE CO_DATA FOR PLOTTING BY MASKING IT
    count_t = daily_data_dict['count_t']
    mask = (count_t == 0) # Masking all zero values in count_t array
    field_mt = ma.array(daily_data_dict['CO_ppb'], mask=mask) # Apply this mask to figtype as well

    
    # 4:
    # # Correct data, for example if:
        # TROPOMI (1) gfed_tropomi (1) & GFED (10) = 12 -> 11
    plumes = daily_data_dict[mask_type]
    plumes[plumes == 12] = 11
    plumes[plumes == 102] = 101
    plumes[plumes == 112] = 111
    
    if subplots:
        SubPlots(lon = lon, lat = lat, data1 = field_mt, data2 = plumes, out_name = out_name)
    elif not subplots:
        CreatePlumesMap(lon, lat, plumes, out_name=out_name_mask)
        CreateCOMap(lon, lat, field_mt, out_name=out_name_co)

    return



def CreatePlumesMap(lon, lat, plumes, out_name):
    """
    Generate png image for plume data ('plumes'), within 'lat' - 'lon' extent, saved at 'out_name'.

    """
    
    plt.figure(figsize=(10,6))#, dpi=1200)
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # Load some Cartopy Features
    ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m')
    ax.add_feature(ocean_50m, edgecolor = 'face', facecolor = '#FFFFFF', zorder=1)#'#d0d0d0', zorder=1) 
    
    # bounds = [no plume, TROPOMI, GFED, TROPOMI+GFED, EDGAR, EDGAR+TROPOMI, GFED+EDGAR, GFED+EDGAR+TROPOMI]
    colors = ['#E5E5E5', '#F77F00', '#8FC93A', '#594236', '#0496FF', '#791E94', '#7AFDD6', '#F0C808']
    cmap = ListedColormap(colors)
    # Setting the (discrete) boundaries of the map
    bounds = [0,1,10,11,100,101,111,112]
    norm = BoundaryNorm(bounds, cmap.N)
    
    ax.pcolormesh(lon, lat, plumes, cmap = cmap, norm=norm, transform=ccrs.PlateCarree())
    it = lambda color: plt.Rectangle((0,0),1,1, facecolor=color, edgecolor='black')
    ax.legend([it(colors[0]), it(colors[1]), it(colors[2]), it(colors[3]), \
                it(colors[4]), it(colors[5]), it(colors[6]), it(colors[7])], \
                       ["no plume", "TROPOMI", "GFED", "TROPOMI + GFED", "EDGAR", \
                        "EDGAR + TROPOMI", "GFED + EDGAR", "GFED + EDGAR + TROPOMI"], \
                           loc='upper center', bbox_to_anchor=(0.5, -0.045), ncol=3, \
                               fancybox=False, shadow=False, frameon=False)
           
    # Save the figure
    plt.savefig(out_name, bbox_inches='tight')#, dpi=1200)
    plt.cla()
    plt.clf()

    return


def CreateCOMap(lon, lat, field_mt, out_name):
    """
    Generate png image for CO data ('field_mt') as measured by TROPOMI, within 'lat' - 'lon' extent, saved at 'out_name'.

    """
    # plot with cartopy
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    

    # Add some cartopy features to the map
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
    ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3) 
        
    cs = plt.pcolormesh(lon, lat, field_mt, cmap='rainbow', transform=ccrs.PlateCarree())
    cbaxes = fig.add_axes([0.2, 0.03, 0.6, 0.03]) 
    cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
    cb.set_label('ppb')
    
    # Save the figure
    fig.savefig(out_name, bbox_inches='tight')#, dpi=1200)
    plt.cla()
    plt.clf()
    #plt.close()
    return


  
def SubPlots(lon, lat, data1, data2, out_name):
    """
    Creates two subplots (CO total column as observed by TROPOMI) and one of the plumes identified in it, compared to the GFED and EDGAR databases.
    """

    plt.ioff()
    
    fig = plt.figure(figsize=(20,6))#, dpi=1200)
    
    # ==============================
    # Plot the CO_ppb colormap
    # ==============================
    ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
    #ax1.coastlines()
    gl = ax1.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # Add some cartopy features to the map
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
    ax1.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3) 
    
    cs = plt.pcolormesh(lon, lat, data1, cmap='rainbow', transform=ccrs.PlateCarree())
    cbaxes = fig.add_axes([0.15, -0.01, 0.3, 0.03]) 
    cb1 = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal')#, fraction=0.046, pad=0.04) 
    cb1.set_label('ppb')
    
    # ==============================
    # Plot the plume mask discrete color map
    # ==============================
    ax2 = plt.subplot(122, projection=ccrs.PlateCarree())
    #ax2.coastlines()
    gl = ax2.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m')
    ax2.add_feature(ocean_50m, edgecolor = 'face', facecolor = '#FFFFFF', zorder=1)#'#d0d0d0', zorder=1) 
    
    # bounds = [no plume, TROPOMI, GFED, TROPOMI+GFED, EDGAR, EDGAR+TROPOMI, GFED+EDGAR, GFED+EDGAR+TROPOMI]
    #colors = ['#16060C', '#FF7621', '#FEB504', '#9F5244', '#083554', '#70CED0', '#FFFFFF', '#FFFFFF']
    #colors = ['#f2f2f2', '#FF7621', '#33a02c', '#663300', '#1f78b4', '#660099', '#FFFFFF', '#FFFFFF']#'‎#0D98BA', '#FFFFFF']
    colors = ['#E5E5E5', '#F77F00', '#8FC93A', '#594236', '#0496FF', '#791E94', '#7AFDD6', '#F0C808']
    cmap = ListedColormap(colors)
    # Setting the (discrete) boundaries of the map
    bounds = [0,1,10,11,100,101,111,112]
    norm = BoundaryNorm(bounds, cmap.N)
    
    ax2.pcolormesh(lon, lat, data2, cmap = cmap, norm=norm, transform=ccrs.PlateCarree())
    it = lambda color: plt.Rectangle((0,0),1,1, facecolor=color, edgecolor='black')
    ax2.legend([it(colors[0]), it(colors[1]), it(colors[2]), it(colors[3]), \
                it(colors[4]), it(colors[5]), it(colors[6]), it(colors[7])], \
                       ["no plume", "TROPOMI", "GFED", "TROPOMI + GFED", "EDGAR", \
                        "EDGAR + TROPOMI", "GFED + EDGAR", "GFED + EDGAR + TROPOMI"], \
                           loc='upper center', bbox_to_anchor=(0.5, -0.045), ncol=3, \
                               fancybox=False, shadow=False, frameon=False)
           
    # Save the figure
    plt.savefig(out_name, bbox_inches='tight')#, dpi=1200)
    plt.close()
    
    return    



def CreateWindVector(daily_data_dict, basepath, skip=30):
        
    """
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
        Contains at least: lat_min, lat_max, lon_min, lon_max, day, month, year,
        and np.array with values to plot and count_t.
    basepath : str
        String to folder where output files will be stored.
    skip : int, optional
        Every 'skip' arrow is drawn in figure. The default is 30.

    Returns
    -------
    Figure(s) stored as png in '...basepath...//04_output//plume_figures//'

    """
    
    # Creating output directory
    out_dir = ut.DefineAndCreateDirectory(os.path.join(basepath + os.path.sep + r'04_output\plume_figures'))
    
    # Creating output filename
    curr_time = ut.GetCurrentTime()
    tropomi_datestamp = f"{daily_data_dict['month']}_{daily_data_dict['day']}_{daily_data_dict['year']}"
    current_timestamp = f"{curr_time['year']}{curr_time['month']}{curr_time['day']}{curr_time['hour']}{curr_time['minute']}{curr_time['second']}"
    out_name = out_dir + rf"fig_subplots_{tropomi_datestamp}___{current_timestamp}.png"
    
    # Define U and V wind 
    U = daily_data_dict['u_wind']
    V = daily_data_dict['v_wind']
    
    #windspeed = (U ** 2 + V ** 2) ** 0.5
    #direction = 180 + (180/np.pi)*np.arctan2(V,U)
    
    # Create a longitude and latitude np.meshgrid in target resolution
    nlat_t, nlon_t = daily_data_dict['u_wind'].shape

    # Creating lat/lon meshgrid, to correctly place grid cells
    lon_t = np.linspace(daily_data_dict['lon_min'], daily_data_dict['lon_max'], nlon_t)
    lat_t = np.linspace(daily_data_dict['lat_min'], daily_data_dict['lat_max'], nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    # Prepare data for plotting
    count_t = daily_data_dict['count_t']
    mask = (count_t == 0) # Masking all zero values in count_t array
    field_mt = ma.array(daily_data_dict['CO_ppb'], mask=mask) # Apply this mask to figtype as well

    # plot with cartopy
    plt.ioff()
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # Add some cartopy features to the map
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m') 
    ax.add_feature(land_50m, edgecolor='k',linewidth=0.5,facecolor='None',zorder=3) 

    # Draw vectors
    cs = plt.quiver(lon[::skip, ::skip], lat[::skip, ::skip], U[::skip, ::skip], V[::skip, ::skip], transform=ccrs.PlateCarree(), pivot = 'mid', zorder=2)
    
    # Use total column CO as background of wind vectors
    cs = plt.pcolormesh(lon, lat, field_mt, cmap='rainbow', transform=ccrs.PlateCarree(), zorder = 1)
    cbaxes = fig.add_axes([0.2, 0.03, 0.6, 0.03]) 
    cb = plt.colorbar(cs, cax = cbaxes, orientation = 'horizontal' )
    cb.set_label('CO concentration (ppb)')
    
    # Save the figure
    plt.savefig(out_name, bbox_inches='tight')#, dpi=1200)
    plt.close()

    return
   
        
        
        
