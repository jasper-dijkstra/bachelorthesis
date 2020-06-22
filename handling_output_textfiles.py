# -*- coding: utf-8 -*-
"""
Created on Sat May 30 09:34:30 2020

@author: Jasper Dijkstra


"""

import os
import numpy as np
from scipy import ndimage

# Local imports
import utilities as ut


def GetStats(daily_data, labeled_plumes, num_features):
    
    # Make a copy of the identified plumes
    plumes = np.copy(daily_data['plumes_explained']) # Get the plumes
    
    # Initiate list with plume types
    plumetypes = []
    
    for i in range(1,num_features+1): # Loop over all labelled plumes
        plume = np.copy(labeled_plumes)
        plume[plume != i] = 0 # Keep only the labelled plume
        
        # Retrieve plume value
        plume_y, plume_x = np.where(plume != 0)
        plumevalue = plumes[plume_y[0], plume_x[0]]
        
        plumetypes.append(int(plumevalue))

    # Get statistics from this information
    daily_data['explained plumes'] = len(np.array(plumetypes)[np.array(plumetypes) != 1])
    daily_data['explained by gfed'] = len(np.array(plumetypes)[np.array(plumetypes) == 11])
    daily_data['explained by edgar'] = len(np.array(plumetypes)[np.array(plumetypes) == 101])
    
    return daily_data


def NotePlumeCoordinates(daily_data_dict, basepath, lonres, latres, params, compare_gfed_edgar):
    """
    
    Parameters
    ----------
    daily_data_dict : dictionary
        daily_data[<day>], contains data about TROPOMI measurement per day.
    TODO!!!!!!!!!!!!!!ct (txt file) of plume coordinates will be stored.

    Returns
    -------
    Saves .txt file in coord_directory.

    """
    
    # Output directories, for coordinates and figures
    coord_directory = ut.DefineAndCreateDirectory(os.path.join(basepath + r'\04_output\plume_coordinates'))
    
    
    # Defining boundaries
    lat_min = daily_data_dict['lat_min']
    lat_max = daily_data_dict['lat_max']
    lon_min = daily_data_dict['lon_min']
    lon_max = daily_data_dict['lon_max']
    
    # Deciding on the nlon_t and nlat_t
    field_t = daily_data_dict['CO_ppb']
    nlon_t = len(field_t[0])
    nlat_t = len(field_t)
    
    if compare_gfed_edgar:
        # Define the plumes
        plumes = daily_data_dict['plumes_explained']
        
        # Label the plumes
        labels, nlabels = ndimage.label(plumes)
        
        # Get some statistics
        daily_data_dict = GetStats(daily_data_dict, labels, nlabels)
    
    else:
        # Define the plumes
        plumes = daily_data_dict['plume_mask']
        
        # Label the plumes
        labels, nlabels = ndimage.label(plumes)
    

    # Get some label statistics, for file
    max_xco = ndimage.measurements.maximum_position(field_t, labels, np.arange(nlabels) + 1)
    mean_xco_raw = ndimage.measurements.mean(field_t, labels, np.arange(nlabels) + 1)
    mean_xco = [round(num, 3) for num in mean_xco_raw]
    plume_size = ndimage.labeled_comprehension(plumes, labels, np.arange(nlabels) + 1, np.sum, float, 0)
    unexplained = round((100 - (daily_data_dict['explained plumes'] / nlabels)*100), 2)
    exp_by_gfed = round(((daily_data_dict['explained by gfed'] / nlabels)*100), 2)
    exp_by_edgar = round(((daily_data_dict['explained by edgar'] / nlabels)*100), 2)

    
    # Generate coordinate meshgrid
    lon_t = np.linspace(lon_min, lon_max, nlon_t)
    lat_t = np.linspace(lat_min, lat_max, nlat_t)
    lon, lat = np.meshgrid(lon_t, lat_t)
    
    # Write to txt file
    day = daily_data_dict['day']
    month = daily_data_dict['month']
    year = daily_data_dict['year']
    curr_time = ut.GetCurrentTime()
    total = daily_data_dict['total_plumes']
    bufferarea = round((np.pi*((params[0]*((lonres+latres)/2))**2)), 1)
    filename = os.path.join(coord_directory + \
        'Plume_coordinates_{}_{}_{}.txt'.format(month, day, year))
    
    headerstring = f"""#----------------------------------------
#----------------------------------------
This file was automatically generated at: {curr_time['year']}/{curr_time['month']}/{curr_time['day']} {curr_time['hour']}:{curr_time['minute']}

This file contains a list with information on Carbon Monoxide plumes at {month}/{day}/{year}, between:
longitudes: [{lon_min}, {lon_max}] 
latitudes: [{lat_min}, {lat_max}] 

column descriptions:
- type:         Plume centre of mass origin: (Unknown, TROPOMI, TROPOMI+GFED, TROPOMI+EDGAR, TROPOMI+GFED+EDGAR)
- latitude:     Latitude of northernmost edge of plume
- longitude:    Longitude of westermost edge of plume
- grid_cells:   Amount of grid cells (~{lonres}x{latres}km) in plume
- CO_max:       Highest Carbon Monoxide concentration measured in plume (ppb)
- CO_average:   Average Carbon Monoxide concentration measured in plume (ppb)


Total amount of plumes identified by TROPOMI: {total}

Percentage of TROPOMI plumes that can be explained by GFED: {exp_by_gfed}%
Percentage of TROPOMI plumes that can be explained by EDGAR: {exp_by_edgar}%
Percentage of TROPOMI plumes that cannot be explained by EDGAR or GFED:{unexplained}%

Note: a TROPOMI grid cell can be explained by EDGAR or GFED when an EDGAR or GFED enhancement was detected within
a circular buffer with radius {params[0]} (~{bufferarea} km^2) around a TROPOMI plume.

Other parameter settings:
    Moving window size:             {params[2]}x{params[2]} grid cells
    Moving window step size:        {params[3]} grid cells
    Standard deviations treshold:   {params[1]}

#----------------------------------------
#----------------------------------------
plume_origin; latitude; longitude; grid_cells; CO_max; CO_average;

"""
    
    f = open(filename, 'w+')
    f.write(headerstring)

    for i in range(1,nlabels+1):
        plume = np.copy(labels)
        plume[plume != i] = 0 # Keep only the labelled plume
        
        # Retrieve plume value
        y, x = np.where(plume != 0) # Get the indices of the left top corner of plume
        y = y[0] # Get only northernmost grid cell coordinate
        x = x[0] # Get only westermost grid cell coordinate
        
        x_co_max = int(max_xco[i-1][1])
        y_co_max = int(max_xco[i-1][0])
        plume_origin = 'Unknown' if (plumes[y,x] == 1) else 'Wildfire' \
            if (plumes[y,x] == 11) else 'Anthropogenic' if (plumes[y,x] == 101) \
                else 'Wildfire/anthropogenic' if (plumes[y,x] == 111) else 'Error'
        f.write(f"{plume_origin}; {lat[y, x]}; {lon[y, x]}; {plume_size[i-1]}; {round(field_t[y_co_max, x_co_max], 3)}; {mean_xco[i-1]}\n")
    f.close()
    
    # Notify a text file has been generated
    #print(f'Generated: {filename}')
    
    return