# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 14:23:54 2020

@author: jaspd

This script:
    1. creates monthly file of data within specified location using (splitted) csv files
    2. splits this monthly file into individual daily csv files
    
"""

import os
from datetime import datetime
import handling_input as inpt


# Set the Target Boundaries (degrees)
lon_min = 100
lon_max = 160
lat_min = -50
lat_max = 0
boundaries = [lat_min, lat_max, lon_min, lon_max]


#-----------------------------------
# 1. MERGING SPLITTED CSV FILES TO ONE, FILTERING DATA THAT FALLS WITHIN EXTENT
#-----------------------------------

print('started merging csv files at: {}'.format(datetime.now()))
start = datetime.now()

# Creating list with all csv files to append
basepath = r'C:/Users/jaspd/Desktop/THESIS_WORKINGDIR/November_2018_RAW/'
csvs = [os.path.join(basepath + 'xaa.csv'), os.path.join(basepath + 'xab.csv'), os.path.join(basepath + 'xac.csv'), os.path.join(basepath + 'xad.csv'), os.path.join(basepath + 'xae.csv'), os.path.join(basepath + 'xaf.csv')]

# Out CSV path
out_csv_file = os.path.join(basepath + 'tropomi_co_Australia_nov18.csv')

# Filter all data that falls within bbox
filtered = []
for file in csvs:
    out_data = inpt.filter_csv_by_bbox(file, boundaries)
    filtered.append(out_data)

# Create one list that contains all grid cell-lists
merged_lists = []
for sublist in filtered:
    for cell in sublist:
        merged_lists.append(cell)

inpt.export_as_csv(out_csv_file, merged_lists)
print('finished merging csv files at: {}'.format(datetime.now()))


#-----------------------------------
# 2. SEPARATING FILTERED DATA TO INDIVIDUAL FILES PER DAY
#-----------------------------------

print('started filtering daily data at: {}'.format(datetime.now()))

# Creating list with all days of the month
days = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]


for day in days:
    csv_day = inpt.filter_csv_by_day(out_csv_file, day)
    exp_path = os.path.join(basepath + 'tropomi_co_Australia_{}_nov18.csv'.format(day))
    inpt.export_as_csv(exp_path, csv_day)

print('finished filtering daily data at: {}'.format(datetime.now()))
print('total time elapsed: {}'.format(datetime.now()-start))