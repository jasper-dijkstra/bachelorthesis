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


def GetStats(daily_data_dict):
    plumes = np.copy(daily_data_dict['plume_mask'].flatten())
    plumes = plumes[plumes > 0]
    frequency = np.bincount(plumes)
    
    # Make sure frequency has got enough indices to complete this function:
    append_values = 114 - len(frequency)
    frequency = np.lib.pad(frequency, ((0,append_values)), 'constant', constant_values=(0))

    total_tropomi = frequency[1] + frequency[12] + frequency[112] + frequency[113] 
    
    plumes[plumes == 12] = 11
    plumes[plumes == 102] = 101
    plumes[plumes == 112] = 111
# =============================================================================
#     plumes = daily_data_dict['plume_mask'].flatten()
#     plumes = plumes[plumes > 0]
#     frequency = np.bincount(plumes)
# =============================================================================

    total_tropomi_in_buffer = frequency[1] + frequency[11] + frequency[111]
    total_tropomi_in_gfed = frequency[11]
    total_tropomi_in_edgar = frequency[101]
    total_tropomi_in_edgar_gfed = frequency[111]

    
    exp_by_gfed = round(100 * ((total_tropomi_in_gfed + total_tropomi_in_edgar_gfed) / total_tropomi_in_buffer), 1)
    exp_by_edgar = round(100 * ((total_tropomi_in_edgar + total_tropomi_in_edgar_gfed) / total_tropomi_in_buffer), 1)
    unexplained = round(100 - exp_by_gfed - exp_by_edgar, 1)
    
    stats = [total_tropomi, exp_by_gfed, exp_by_edgar, unexplained]

    return stats


def NotePlumeCoordinates(daily_data_dict, basepath, lonres, latres):
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
    
    # Getting statistics on the plumes
    stats = GetStats(daily_data_dict)
    
    # Deciding on the nlon_t and nlat_t
    field_t = daily_data_dict['CO_ppb']
    nlon_t = len(field_t[0])
    nlat_t = len(field_t)
    
    # Group plumes, so they will not be identified as independent cells
    plumes = daily_data_dict['plume_mask']
    
    # Correct the data:
        # Only TROPOMI plumes must be identified, therefore:
    plumes[plumes == 10] = 0    # Remove GFED
    plumes[plumes == 11] = 0    # Remove GFED within TROPOMI buffer
    plumes[plumes == 100] = 0   # Remove EDGAR
    plumes[plumes == 101] = 0   # Remove EDGAR within TROPOMI buffer
    plumes[plumes == 110] = 0   # Remove EDGAR + GFED overlap
    # Now only TROPOMI plumes whose grid cells directly overlap with EDGAR or GFED-
    # are left. GFED and EDGAR within TROPOMI bufferzones are deleted
    # Values: [0 = No plume, 1 = TROPOMI, 12 = TROPOMI & GFED, 102 = TROPOMI & EDGAR
    # 112 = TROPOMI & GFED & EDGAR]
    
    # Now label the data
    plumes_0 = plumes > 0 # Define data to be labelled
    labels, nlabels = ndimage.label(plumes_0) # Label data
    indices = ndimage.center_of_mass(plumes, labels, np.arange(nlabels) + 1) # Define center of mass of labeled plumes
    
    # Get some label statistics
    max_xco = ndimage.measurements.maximum_position(field_t, labels, np.arange(nlabels) + 1)
    mean_xco_raw = ndimage.measurements.mean(field_t, labels, np.arange(nlabels) + 1)
    mean_xco = [round(num, 3) for num in mean_xco_raw]
    plume_size = ndimage.labeled_comprehension(plumes, labels, np.arange(nlabels) + 1, np.sum, float, 0)
    
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
- latitude:     Latitude of center of plume
- longitude:    Longitude of center of plume
- grid_cells:   Amount of grid cells (~{lonres}x{latres}km) in plume
- CO_max:       Highest Carbon Monoxide concentration measured in plume (ppb)
- CO_average:   Average Carbon Monoxide concentration measured in plume (ppb)


Total amount of plumes identified by TROPOMI: {total}

Percentage of TROPOMI plumes that can be explained by GFED: {stats[1]}%
Percentage of TROPOMI plumes that can be explained by EDGAR: {stats[2]}%
Percentage of TROPOMI plumes that cannot be explained by EDGAR or GFED:{stats[3]}%

Note: a TROPOMI grid cell can be explained by EDGAR or GFED when an EDGAR or GFED enhancement was detected within
a circular buffer with radius 7 (~70 km) around a TROPOMI plume.

#----------------------------------------
#----------------------------------------
type; latitude; longitude; grid_cells; CO_max; CO_average;

"""
    
    f = open(filename, 'w+')
    f.write(headerstring)

    for i in range(len(indices)):
        x = int(indices[i][1])
        y = int(indices[i][0])
        x_co_max = int(max_xco[i][1])
        y_co_max = int(max_xco[i][0])
        plume_origin = 'TROPOMI' if (plumes[y,x] == 1) else 'TROPOMI+GFED' \
            if (plumes[y,x] == 12) else 'TROPOMI+EDGAR' if (plumes[y,x] == 102) \
                else 'TROPOMI+GFED+EDGAR' if (plumes[y,x] == 112) else 'Unknown origin'
        f.write(f"{plume_origin}; {lat[y, x]}; {lon[y, x]}; {plume_size[i]}; {round(field_t[y_co_max, x_co_max], 3)}; {mean_xco[i]}\n")
    f.close()
    
    print(f'Generated: {filename}')
    
    return