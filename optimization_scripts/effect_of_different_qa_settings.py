# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 19:18:14 2020

@author: Jasper Dijkstra

Script that calculates the effect of applying different qa-filter standards:
    how many grid cells are not considered when applying filter x?

standards:
    ESA-standard: qa > 0.5
    model-standard: 0.1 < qa < 2.0

"""

import pandas as pd
import numpy as np
import os, sys


# Local imports
sys.path.append("..") # get one directory up
import utilities as ut


inputDirectory = r'C:\Users\jaspd\Documents\Python\00_bachelorthesis\bachelorthesis\THESIS_WORKINGDIR\00_daily_csv'
files = ut.ListFilesInDirectory(inputDirectory, extension='.csv', maxfiles=None)

bbox = [-50, 0, 100, 160]

# Dataframe header values
CO_header = ['lat0', 'lat1', 'lat2', 'lat3', 'lat', 'lon0', 'lon1', 'lon2', 'lon3', 'lon', 'xco', 'xco_unc', 'xco_ppb', 'xco_ppb_unc', 'qa', 'weekday', 'day', 'month', 'aerosol_opthick', 'aerosol_layer', 'orbitnr', 'time']

esa_standard = []
model_standard = []
for filtertype in [1, 2]:
    for csv_file in files:
        # Read the csv file as Pandas DataFrame
        COdata = pd.read_csv(csv_file, header=None, names=CO_header)
        
        # Identify if observations falls within latlon range
        COdata['valid'] = (COdata['lat'] >= bbox[0]) & (COdata['lat'] < bbox[1]) \
            & (COdata['lon'] >= bbox[2]) & (COdata['lon'] < bbox[3])
        
        # Remove all values that do not fall within latlon range
        COdata = COdata.drop(COdata[COdata['valid'] == False].index)
        del COdata['valid']
        
        # Converting data to numpy array
        COdata = COdata.to_numpy() # Transforming pd.dataframe to np.array
        xco_ppb_1d = COdata[:,12] # Get CO data from array
        qa_1d = COdata[:,14] # Get qa data from array
        
        # Give indices with a too high a 'nodata' value
        if filtertype == 1:
            qa_indices = np.where((qa_1d <= 0.1) | (qa_1d >= 2.0))[0]
        if filtertype == 2:
            qa_indices = np.where((qa_1d < 0.5))[0]
        xco_ppb_1d[qa_indices] = np.nan
        
        # What percentage of data is removed when applying filter?
        percentage_removed = (1 - np.count_nonzero(~np.isnan(xco_ppb_1d)) / len(xco_ppb_1d))

        # Append this to list
        if filtertype == 1:
            model_standard.append(percentage_removed)
        if filtertype == 2:
            esa_standard.append(percentage_removed)

# Mean deviation (with error margin (one st.dev)) per day from model
mean_model = np.mean(np.array(model_standard))
mean_model_margin = np.nanstd(np.array(model_standard))

# Mean deviation (with error margin (one st.dev)) per day from ESA
mean_esa = np.mean(np.array(esa_standard))
mean_esa_margin = np.nanstd(np.array(esa_standard))
