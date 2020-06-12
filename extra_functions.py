# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:16:24 2020

@author: Jasper Dijkstra

Script contains function that are no longer used in TROPOMI plume detection algorithm:
    - CheckCO() - Check CO concentration (in ppb) for given lat/lon combination
    - SetMostOccuring() - Define most occuring value in array, and assign this value to all cells

"""

import numpy as np
import collections


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


def SetMostOccuring(arr):
    """
    Define most occuring value in array, and assign this value to all cells

    If some values occur just as much as others, the median value will be assigned

    """
    
    one_dimension = arr.flatten()
    counts = collections.Counter(one_dimension).most_common(2)
    #bincount = np.bincount(one_dimension)
    #idx = np.where(bincount == np.max(bincount))[0]
    
    # Make sure 0 is not most occuring
    try:
        if counts[0][0] == 0:
            one_dimension[:] = counts[1][0]
        else:
            one_dimension[:] = counts[0][0]
    except IndexError:
        one_dimension[:] = counts[0][0]

    
    outarr = one_dimension.reshape(arr.shape)
    
    return outarr

