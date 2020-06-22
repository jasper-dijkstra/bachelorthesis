# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:28:24 2020

@author: Jasper Dijkstra

This algorithm detects Carbon Monoxide plumes in TROPOMI data

Edit the parameters section to tune the algorithm as desired.

The rest of the script is built up as follows:
    - INITIALIZATION:
        TROPOMI data is processed to np.arrays
        If desired, measurements above oceans can be filtered out 
        Data is stored in dictionaries (daily_data_dict)
    - MAIN MODEL:
        Using a moving window, highest values are separated from background
        TROPOMI plumes can be compared to data from the Global Fire Emissions Database
    - HANDLING OUTPUT:
        If desired, textfiles containing plume coordinates are generated
        If desired, figures are generated

"""

#--------------------
# IMPORTING FUNCTIONS
#--------------------
import os
import copy
from datetime import datetime
from scipy import ndimage
import numpy as np

# Local imports
import handling_input as inpt
import handling_GFED as GFED
import handling_EDGAR as EDGAR
import handling_output_textfiles as txt
import handling_output_figures as figs
import masking_functions as mask
import raster_tools as raster
import fetching_CAMS as cams
import fetching_winddata as wind
import utilities as ut


# ========================================================
# NON-USER DEFINED PARAMETERS
# ========================================================

def InitializingParameters(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath):
    # Store all defined extents in a list:
    boundaries = [lat_min, lat_max, lon_min, lon_max]
    
    # Calculating resolution based on lonres, latres
    target_lon = int(abs((lon_max-lon_min)/(lonres/110)))
    target_lat = int(abs((lat_max-lat_min)/(latres/110)))
    
    # Create a list with all files to apply the analysis on
    input_files_directory = os.path.join(basepath + r'00_daily_csv\\')
    files = ut.ListFilesInDirectory(input_files_directory, maxfiles=4)
    del files[0:3] # Make sure script only runs for October 13, 2018 #DELETE THIS LATER

    return boundaries, target_lon, target_lat, files


# ========================================================
# Initialization (reading data)
# ========================================================

def ValidateInputs(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath, \
                   GFED_path, EDGAR_path, CAMS_path, behaviour_settings):
    """
    Assert user defined parameters will not raise errors later in the algorithm
    """
    # Assert sure extents fall within boundary
    assert -180 <= lon_min <= 180 or -180 <= lon_max <= 180 or lon_min <= lon_max, 'maximum longitude cannot be smaller than or equal to minimum, and should be within range (-180, 180)!'
    assert -90 <= lat_min <= 90 or -90 <= lat_max <= 90 or lat_min <= lat_max, 'maximum latitude cannot be smaller than or equal to minimum, and should be within range (-90, 90)!'
    
    # Assert resolution is larger than TROPOMI minimum:
    assert lonres > 7, 'TROPOMI minimum longitude resolution is 7 km!'
    assert latres > 7, 'TROPOMI minimum latitude resolution is 7 km!'
    
    # Assert if given directories exist
    if behaviour_settings[1] == True:
        assert os.path.isdir(CAMS_path), f'Directory {CAMS_path} was not found!'
    if behaviour_settings[2] == True:
        assert os.path.isdir(GFED_path), f'Directory {GFED_path} was not found!'
        assert os.path.isdir(EDGAR_path), f'Directory {EDGAR_path} was not found!'

    
    return


def CollectingData(boundaries, target_lon, target_lat, files, basepath, \
                   CAMS_path, apply_land_sea_mask, use_wind_rotations, \
                       incorporate_cams):
    """
    Opens the data required to execute plume detection algorithm

    Parameters
    ----------
    boundaries : list, tuple
        list, tuple with desired data extents [lat_min, lat_max, lon_min, lon_max]
    target_lon : int, target latitude (grid cells)
    target_lat : int, target longitude (grid cells)
    files : list
        List with all csv files containing TROPOMI data in directory.
    basepath : str
        Path to working directory (outputs will be stored here).
    CAMS_path : str
        Path to folder where CAMS files are stored
    apply_land_sea_mask : bool
    use_wind_rotations : bool
    incorporate_cams : bool

    Returns
    -------
    daily_data : dict
        Dictionary with all data required to execute algorithm, sorted per day.

    """
    # Setting the time of starting the script
    start = datetime.now()
    
    # Reading daily csv files for specified area and day as np.arrays
    daily_data = {}
    for i, file in enumerate(files):    
        # Reading daily csv's as input array
        daily_data[i] = inpt.CSVtoArray(file, boundaries, target_lon, target_lat)
        
        # Remove background, by CAMS observations
        if incorporate_cams:
            dates = [daily_data[i]['day'], daily_data[i]['month'], daily_data[i]['year']]
            bbox = [daily_data[i]['lat_min'], daily_data[i]['lat_max'], daily_data[i]['lon_min'], daily_data[i]['lon_max']]
            xres = int((110 * (bbox[3]-bbox[2])) / len(daily_data[i]['CO_ppb'][0]))
            yres = int((110 * (bbox[1]-bbox[0])) / len(daily_data[i]['CO_ppb']))
            cams_arr = cams.FetchCams(CAMS_path, dates, bbox, xres, yres)
            daily_data[i]['CO_excl_background'] = daily_data[i]['CO_ppb'] - cams_arr
        
        # Filter measurements taken above the oceans (higher uncertainty)
        if apply_land_sea_mask:
            daily_data[i]['CO_ppb'] = mask.land_sea_mask(daily_data[i]['CO_ppb'], boundaries)
            daily_data[i]['count_t'] = mask.land_sea_mask(daily_data[i]['count_t'], boundaries)
    
        # collect meteodata via ECMWF CDS API:
        if use_wind_rotations:
            u_wind, v_wind = wind.FetchWindData(daily_data[i], pressure=700, timerange=6, basepath=basepath)
            daily_data[i]['u_wind'] = u_wind
            daily_data[i]['v_wind'] = v_wind
    
    print('Total time elapsed reading data: {}'.format(datetime.now()-start))

    return daily_data


# ========================================================
# DETECTION ALGORITHM
# ========================================================

def Detection(params, daily_data, boundaries, GFED_path, EDGAR_path, lonres, latres, \
              apply_overlap_filter, use_wind_rotations, incorporate_cams, compare_gfed_edgar):
    
    # Note starting time of function
    start_main = datetime.now()
    
    """ Check for plumes based on four criteria: """
    for day in daily_data:
        """ 1: At least 2 standard deviations above average of moving window """
        # Apply moving Window operation, with copy of CO_ppb
        if incorporate_cams:
            arr = copy.deepcopy(daily_data[day]['CO_excl_background'])
        else:
            arr = copy.deepcopy(daily_data[day]['CO_ppb'])
        
        plumes, co_average = raster.MovingWindow(arr, mask.identify_enhancements, \
            window = (params[2],params[2]), step = params[3], st_devs = params[1]) 
        
        # Check if plume was detected in at least 95% of windows
        plumes[plumes >= 0.95] = 1 # If true, plume (1)
        plumes[plumes < 0.95] = 0 # If false, no plume (0)
        
        # Removing nuisances by making sure there is at least one neighbour
        neighbors = raster.CountNeighbors(plumes)  # Identify neighbors of each grid cell
        plumes[neighbors <= 1] = 0 # If there are 1 or fewer neighbors, undo identification as plume
        
        
        """ 2: Check if plumes overlap > 3 days """
        if apply_overlap_filter:
            print('Overlap filter not available!')


        """ 3: Check if plumes overlap with modelled GFED and EDGAR data, by drawing a buffer around TROPOMI data """
        if compare_gfed_edgar:
            plumes_output = GFEDEDGARComparison(plumes, daily_data[day], boundaries, lonres, latres,\
                                         GFED_path, EDGAR_path, buffersize = params[0])
            daily_data[day]['plume_mask'] = plumes_output[0]
            daily_data[day]['plumes_explained'] = plumes_output[1]
            daily_data[day]['plumes_incl_gfed_edgar'] = plumes_output[2]
        else:
            # Append to daily data dictionary
            daily_data[day]['plume_mask'] = plumes
        
        """ 4: Rotate in the wind, to check if plume is detected downwind and not upwind """
        if use_wind_rotations:
            print('Wind rotations not available!')
            # Rotate identified plumes so they will be in line with the wind direction
            # If plume occurs exclusively downwind, it can be validated as plume.

        print('Total time elapsed executing plume masking algorithm: {}'.format(datetime.now()-start_main))

    return daily_data


# =============================================================================
# def GFEDEDGARComparison2(plumes, daily_data, boundaries, lonres, latres, \
#                         GFED_path, EDGAR_path, buffersize):
#     # Read GFED and EDGAR data
#     gfed_array = GFED.OpenGFED(GFED_path, boundaries, daily_data['day'], \
#                                daily_data['month'], daily_data['year'], lonres, latres)
#     gfed_array[gfed_array > 0] = 1 # As we are only interested to know if emissions happened, return emissions true(1)/false(0)
#     
#     edgar_file = 'v432_CO_2012.0.1x0.1.nc' # EDGAR indsutry filename 
#     edgar_array = EDGAR.OpenEDGAR(os.path.join(EDGAR_path + edgar_file), boundaries, lonres, latres)
#     edgar_array[edgar_array > 0] = 1 # As we are only interested to know if emissions happened, return emissions true(1)/false(0)
# 
#     
#     gfed_array[gfed_array > 0] = 10     # GFED plumes = 10
#     edgar_array[edgar_array > 0] = 100  # EDGAR plumes = 100
#     
#     # Draw Buffer around TROPOMI plumes
#     plumes_buffered = raster.DrawCircularBuffer(plumes, radius = buffersize) # Draw circular buffers around TROPOMI plumes
# 
#     # Delete GFED and EDGAR values outside of plume buffer
#     gfed_tropomi = (plumes_buffered * gfed_array) # GFED within TROPOMI
#     gfed_tropomi[gfed_tropomi > 0] = 1
#     edgar_tropomi = (plumes_buffered * edgar_array) # EDGAR within TROPOMI
#     edgar_tropomi[edgar_tropomi > 0] = 1
#     
#     # Adding all arrays with different plume origins:
#     plumes = (plumes + gfed_array + edgar_array + gfed_tropomi + edgar_tropomi).astype(int)
#     
#     # Now: (0: no plume, 1: TROPOMI plume, 10: GFED plume (within bufferzone of TROPOMI plume, \
#        # 11: TROPOMI + GFED identified plume))
# 
#     return plumes
# =============================================================================


def GFEDEDGARComparison(plumes, daily_data, boundaries, lonres, latres, \
                        GFED_path, EDGAR_path, buffersize):
    
    # Read GFED and EDGAR data
    gfed_array = GFED.OpenGFED(GFED_path, boundaries, daily_data['day'], \
                               daily_data['month'], daily_data['year'], lonres, latres)
    gfed_array[gfed_array > 0] = 1 # As we are only interested to know if emissions happened, return emissions true(1)/false(0)
    
    edgar_file = 'v432_CO_2012.0.1x0.1.nc' # EDGAR indsutry filename 
    edgar_array = EDGAR.OpenEDGAR(os.path.join(EDGAR_path + edgar_file), boundaries, lonres, latres)
    edgar_array[edgar_array > 0] = 1 # As we are only interested to know if emissions happened, return emissions true(1)/false(0)

    gfed_array[gfed_array > 0] = 10     # GFED plumes = 10
    edgar_array[edgar_array > 0] = 100  # EDGAR plumes = 100
    
    # Initiate array with all 'explained' plumes
    identified_plumes = np.zeros(plumes.shape)

    # Label plumes
    labeled_plumes, num_features = ndimage.label(plumes)
    daily_data['total_plumes'] = num_features
    print(f"Total plumes identified for {daily_data['day']}/{daily_data['month']}/{daily_data['year']}: {num_features}")
    
    for i in range(1,num_features+1): # Loop over all labelled plumes
        bufferfield = copy.deepcopy(labeled_plumes)
        bufferfield[bufferfield != i] = 0 # Keep only the labelled plume
        
        # Draw buffer around labelled plume
        buffer = raster.DrawCircularBuffer(arr = bufferfield, radius = buffersize)
        
        # Delete GFED and EDGAR values outside of plume buffer
        gfed_tropomi = (buffer * gfed_array).astype(int) # GFED within TROPOMI
        gfed_tropomi = int(sum(gfed_tropomi.flatten())/10)
        edgar_tropomi = (buffer * edgar_array).astype(int) # EDGAR within TROPOMI
        edgar_tropomi = int(sum(edgar_tropomi.flatten())/100)
        
        # Append GFED and EDGAR explanation value, if GFED or EDGAR falls within buffer
        if gfed_tropomi > 0:
            identified_plumes[bufferfield == i] += 10
        if edgar_tropomi > 0:
            identified_plumes[bufferfield == i] += 100
    
    # Define the three different outputs
    plumes_tropomi = copy.deepcopy(plumes)
    plumes_explained = plumes + identified_plumes
    plumes_incl_gfed_edgar = plumes_explained + gfed_array + edgar_array
    
    # Correct invalid values, caused by directly overlapping values
    plumes_incl_gfed_edgar[plumes_incl_gfed_edgar == 21] = 11 # GFED + TROPOMI
    plumes_incl_gfed_edgar[plumes_incl_gfed_edgar == 201] = 101 # EDGAR + TROPOMI
    plumes_incl_gfed_edgar[plumes_incl_gfed_edgar == 121] = 111 # GFED + EDGAR + TROPOMI
    plumes_incl_gfed_edgar[plumes_incl_gfed_edgar == 211] = 111 # GFED + EDGAR + TROPOMI
    
    # Now: (0: no plume, 1: TROPOMI plume, 10: GFED plume (within bufferzone of TROPOMI plume, \
       # 11: TROPOMI + GFED identified plume))

    return plumes_tropomi, plumes_explained, plumes_incl_gfed_edgar


# ========================================================
# GENERATING OUTPUTS
# ========================================================

def GeneratingOutputs(daily_data, basepath, lonres, latres, params, \
                      use_wind_rotations, gen_fig_wind_vector, compare_gfed_edgar):
    
    start_end = datetime.now()
    
    # Generate desired outputs based on earlier set preferences
    for day in daily_data:
        # Generate text files with plume coordinates
        daily_data_dict = copy.deepcopy(daily_data[day])
        txt.NotePlumeCoordinates(daily_data_dict, basepath, lonres, latres, params, compare_gfed_edgar)
        
        # Generate figures (subplots)
        if compare_gfed_edgar:
            figs.PlotFigures(daily_data[day], basepath, subplots=True)
            figs.PlotFigures(daily_data[day], basepath, subplots=True, mask_type = 'plumes_explained')
            figs.PlotFigures(daily_data[day], basepath, subplots=True, mask_type = 'plumes_incl_gfed_edgar')
        if not compare_gfed_edgar:
            figs.PlotFigures(daily_data[day], basepath, subplots=True)
        
        # Generate figures of the wind vectors
        if gen_fig_wind_vector:
            assert use_wind_rotations == True, 'use_wind_rotations has to be True before wind vectorfield can be created!'
            figs.CreateWindVector(daily_data[day], basepath, skip=30)
    
    
    print('Total time elapsed generating output: {}'.format(datetime.now()-start_end))
    
    return

    
# ========================================================
# EXECUTE EVERYTHING IN CORRECT ORDER
# ========================================================
def PlumeDetection(lat_min, lat_max, lon_min, lon_max, lonres, latres, \
                   basepath, GFED_path, EDGAR_path, CAMS_path, \
                       params, behaviour_settings):
    
    # First check if all inputs are valid
    ValidateInputs(lat_min, lat_max, lon_min, lon_max, lonres, latres, \
                   basepath, GFED_path, EDGAR_path, CAMS_path, behaviour_settings)
    
    # Notify the algorithm has started
    start = datetime.now()
    print(f'Algorithm started at: {start}')
    print('')
    print('====================')
    
    # Redefine model behaviour settings
    apply_land_sea_mask = behaviour_settings[0]
    incorporate_cams = behaviour_settings[1]
    compare_gfed_edgar = behaviour_settings[2]
    apply_overlap_filter = False # TODO!
    use_wind_rotations = False # TODO!
    gen_fig_wind_vector = False
    if use_wind_rotations:
        gen_fig_wind_vector = True
    
    # Initialize some other required parameters
    boundaries, target_lon, target_lat, files = InitializingParameters(lat_min, lat_max, lon_min, lon_max, lonres, latres, basepath)
    
    # Read all required data, and store in dictionary
    daily_data = CollectingData(boundaries, target_lon, target_lat, files, \
                                basepath, CAMS_path, apply_land_sea_mask, \
                                    use_wind_rotations, incorporate_cams)
        
    print('====================')
    print('')
    print('====================')
    
    # Execute algorithm itself
    daily_data = Detection(params, daily_data, boundaries, GFED_path, EDGAR_path, \
                           lonres, latres, apply_overlap_filter, use_wind_rotations, \
                               incorporate_cams, compare_gfed_edgar)
        
    print('====================')
    print('')
    print('====================')
    
    # Generate outputs
    GeneratingOutputs(daily_data, basepath, lonres, latres, params, use_wind_rotations, \
                      gen_fig_wind_vector, compare_gfed_edgar)
    
    print('====================')
    print('')
    print('====================')
    
    # Notify the algorithm is finished, and how much outputs have been generated
    print(f'Algorithm finished at: {start}')
    print(f'Total time elapsed: {datetime.now()-start}')
    print(f'A total of {len(files)} figures and textfiles were generated at: {basepath}')
    print('')
    
    return daily_data
